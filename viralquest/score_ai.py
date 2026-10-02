import json
import re
import time
from enum import Enum
from pathlib import Path

from loguru import logger

from viralquest.biodata import (
    LLM_ERROR_CLASSES,
    NOVELTY_TIERS,
    BlastnResult,
    BlastxResult,
    HmmDomain,
    LlmOutput,
    NucSequence,
    Orf,
)


class LlmMode(Enum):
    # Pfam domain `details` text is never sent to the LLM (removed to avoid
    # token bloat on Pfam-rich sequences); the high/low distinction is in the
    # other dimensions below. In both modes ORFs without any HMM domain are
    # summarised (count + longest lengths) instead of listed one by one.
    HIGH = "high"   # all hits, full family description, annotated ORFs as objects
    LOW  = "low"    # best hits only, standardised family description, ORFs as one-line summaries


# ---------------------------------------------------------------------------
# Serialisation helpers  (sequences always omitted)
# ---------------------------------------------------------------------------

def _blastx_dict(h: BlastxResult) -> dict:
    return {
        "subject":       h.subject_title,
        "species":       h.species,
        "pct_identity":  h.pct_identity,
        "query_coverage": h.query_coverage,
        "aln_length":    h.aln_length,
        "e_value":       h.e_value,
        "bit_score":     h.bit_score,
    }


def _blastn_dict(h: BlastnResult) -> dict:
    return {
        "subject":        h.stitle,
        "pct_identity":   h.pident,
        "query_coverage": h.qcovhsp,
        "e_value":        h.evalue,
    }


# Longest ORFs without any HMM domain reported to the LLM: a long unannotated
# ORF is a real signal for a novel virus; dozens of short ones are not.
_UNANNOTATED_LONGEST_N = 3


def _round_evalue(e: float | None) -> float | None:
    """Two significant digits — 1.2345e-45 → 1.2e-45 (fewer tokens, same meaning)."""
    return None if e is None else float(f"{e:.2g}")


def _domain_dict(d: HmmDomain) -> dict:
    r = {
        "database":    d.database,
        "target":      d.target,
        "score":       round(d.score, 1),
        "e_value":     _round_evalue(d.e_value),
        "description": d.description,
    }
    if d.type:                      # empty for RVDB / Vfam / EggNOG
        r["type"] = d.type
    return r


def _orf_dict(orf: Orf) -> dict:
    return {
        "name":        orf.name,
        "strand":      orf.strand,
        "frame":       orf.frame,
        "length_aa":   orf.length_aa,
        "hmm_domains": [_domain_dict(d) for d in orf.domains],
    }


def _domain_label(d: HmmDomain) -> str:
    """One-line domain for LOW mode, e.g. "Pfam RdRP_1 – RNA-dependent RNA polymerase [152]"."""
    desc = f" – {d.description}" if d.description and d.description != d.target else ""
    return f"{d.database} {d.target}{desc} [{round(d.score)}]"


def _orf_summary_line(orf: Orf, seq_id: str = "") -> str:
    """
    LOW mode: "ORF_10_1240+ (410 aa): Pfam RdRP_1 – … [152]; RVDB FAM0001 – Picornaviridae [98]".
    The sequence-id prefix of the ORF name is dropped (sequence_id is already in the input).
    """
    name = orf.name.removeprefix(f"{seq_id}_") if seq_id else orf.name
    return f"{name} ({orf.length_aa} aa): " + "; ".join(_domain_label(d) for d in orf.domains)


def _unannotated_summary(orfs: list[Orf]) -> dict:
    """ORFs with no HMM domain, summarised instead of listed."""
    lengths = sorted((o.length_aa for o in orfs), reverse=True)
    return {"count": len(lengths), "longest_aa": lengths[:_UNANNOTATED_LONGEST_N]}


def _orfs_payload(orfs: list[Orf], high: bool, seq_id: str = "") -> dict:
    """
    ORF section of the LLM input.

    HIGH — annotated ORFs as full objects (with their HMM domains).
    LOW  — annotated ORFs as one-line summaries with their domain names.
    Both — ORFs without any domain reduced to a count + the longest lengths.
    """
    annotated   = [o for o in orfs if o.domains]
    unannotated = [o for o in orfs if not o.domains]
    if high:
        payload = {"orfs": [_orf_dict(o) for o in annotated]}
    else:
        payload = {"orf_summary": {
            "total":          len(orfs),
            "with_domains":   len(annotated),
            "annotated_orfs": [_orf_summary_line(o, seq_id) for o in annotated],
        }}
    payload["unannotated_orfs"] = _unannotated_summary(unannotated)
    return payload


def _to_user_message(user_json: dict) -> str:
    """Compact JSON (no spaces after separators) — same content, fewer tokens."""
    return json.dumps(user_json, ensure_ascii=False, separators=(",", ":"))


# ---------------------------------------------------------------------------
# PromptBuilder
# ---------------------------------------------------------------------------

class PromptBuilder:
    """Loads system prompts from llm-info/ and builds the per-sequence input JSON."""

    _PROMPT_DIR = Path(__file__).parent / "llm-info"

    @classmethod
    def system_prompt(cls) -> str:
        return (cls._PROMPT_DIR / "system_prompt.json").read_text(encoding="utf-8")

    @classmethod
    def build_input(cls, seq: NucSequence, mode: LlmMode) -> dict:
        high = mode == LlmMode.HIGH

        # --- taxonomy ---
        tax = seq.taxonomy
        taxonomy_dict = None
        if tax:
            if high:
                taxonomy_dict = {
                    "scientific_name": tax.scientific_name,
                    "clade":    tax.clade,
                    "kingdom":  tax.kingdom,
                    "phylum":   tax.phylum,
                    "class":    tax.class_,
                    "order":    tax.order,
                    "family":   tax.family,
                    "genus":    tax.genus,
                    "species":  tax.species,
                    "genome":   tax.genome,
                }
            else:
                taxonomy_dict = {
                    "family":  tax.family,
                    "genus":   tax.genus,
                    "order":   tax.order,
                    "genome":  tax.genome,
                }

        # --- viral family info ---
        vfi = seq.viral_family_info
        family_info = None
        if vfi:
            family_info = {
                "source": vfi.source,
                "type":   vfi.type,
                "name":   vfi.name,
                "info":   vfi.info_high if high else vfi.info_low,
            }

        # --- BLAST hits ---
        blastx_pool = seq.blastx_nr_hits or seq.blastx_hits
        if high:
            blastx_hits = [_blastx_dict(h) for h in blastx_pool]
            blastn_hits = [_blastn_dict(h) for h in seq.blastn_hits]
        else:
            blastx_hits = [_blastx_dict(seq.best_blastx)] if seq.best_blastx else []
            blastn_hits = [_blastn_dict(seq.best_blastn)] if seq.best_blastn else []

        # Salmon TPM / reads — included when available regardless of token mode,
        # because expression level is a direct quality signal for the LLM.
        salmon_quant = None
        if seq.salmon_tpm is not None or seq.salmon_reads is not None:
            salmon_quant = {}
            if seq.salmon_tpm is not None:
                salmon_quant["tpm"]       = round(float(seq.salmon_tpm),   4)
            if seq.salmon_reads is not None:
                salmon_quant["num_reads"] = round(float(seq.salmon_reads), 1)

        return {
            "sequence_id":      seq.id,
            "sequence_stats": {
                "length_nt":   seq.length,
                "gc_content":  seq.gc_content,
                "n_count":     seq.n_count,
            },
            "is_viral_flag":      seq.is_viral,
            "salmon_quantification": salmon_quant,
            "taxonomy":           taxonomy_dict,
            "blastx_hits":        blastx_hits,
            "blastn_hits":        blastn_hits,
            **_orfs_payload(seq.orfs, high, seq.id),
            "viral_family_info":  family_info,
        }


# ---------------------------------------------------------------------------
# LLM backends  (one class per provider)
# ---------------------------------------------------------------------------

class _LlmBackend:
    needs_rate_limit: bool = True

    def call(self, system_prompt: str, user_json: dict) -> str:
        raise NotImplementedError

    def release(self) -> None:
        """Free resources once scoring ends (no-op for cloud APIs)."""


class OllamaBackend(_LlmBackend):
    """Direct Ollama backend. Uses think=False to suppress chain-of-thought."""

    needs_rate_limit = False  # local process, no API rate limits

    def __init__(self, model_name: str):
        self.model_name = model_name

    def release(self) -> None:
        """
        Unload the model right away. Ollama otherwise keeps it in memory for
        5 minutes after the last request — RAM/VRAM the following pipeline
        steps (and the next sample's Diamond/HMMER) need.
        """
        import ollama
        ollama.generate(model=self.model_name, keep_alive=0)
        logger.info(f"Ollama: model '{self.model_name}' unloaded from memory.")

    def call(self, system_prompt: str, user_json: dict) -> str:
        import ollama
        messages = [
            {"role": "system", "content": system_prompt},
            {"role": "user",   "content": _to_user_message(user_json)},
        ]
        try:
            response = ollama.chat(
                model=self.model_name,
                messages=messages,
                think=False,                          # suppress reasoning (ollama ≥0.4.7)
                options={"temperature": 0.3},
            )
        except TypeError:
            # older ollama-python without think parameter
            response = ollama.chat(
                model=self.model_name,
                messages=messages,
                options={"temperature": 0.3},
            )
        return response["message"]["content"]


class OpenAIBackend(_LlmBackend):
    """OpenAI / compatible API backend. Uses json_object response format."""

    def __init__(self, model_name: str, api_key: str):
        self.model_name = model_name
        self.api_key    = api_key

    def call(self, system_prompt: str, user_json: dict) -> str:
        from openai import OpenAI
        client = OpenAI(api_key=self.api_key)
        response = client.chat.completions.create(
            model=self.model_name,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user",   "content": _to_user_message(user_json)},
            ],
            temperature=0.3,
            response_format={"type": "json_object"},  # enforces JSON output
        )
        return response.choices[0].message.content


class AnthropicBackend(_LlmBackend):
    """Anthropic Claude backend. Filters text-only content blocks to skip thinking."""

    def __init__(self, model_name: str, api_key: str):
        self.model_name = model_name
        self.api_key    = api_key

    def call(self, system_prompt: str, user_json: dict) -> str:
        import anthropic
        client   = anthropic.Anthropic(api_key=self.api_key)
        response = client.messages.create(
            model=self.model_name,
            max_tokens=1024,
            temperature=0.3,
            system=system_prompt,
            messages=[{"role": "user", "content": _to_user_message(user_json)}],
        )
        # take only text blocks — skips any extended-thinking blocks automatically
        return "".join(block.text for block in response.content if block.type == "text")


class GoogleBackend(_LlmBackend):
    """Google Gemini backend. Uses json MIME type to enforce structured output."""

    def __init__(self, model_name: str, api_key: str):
        self.model_name = model_name
        self.api_key    = api_key

    def call(self, system_prompt: str, user_json: dict) -> str:
        from google import genai
        from google.genai import types
        client   = genai.Client(api_key=self.api_key)
        response = client.models.generate_content(
            model=self.model_name,
            contents=_to_user_message(user_json),
            config=types.GenerateContentConfig(
                system_instruction=system_prompt,
                temperature=0.3,
                response_mime_type="application/json",
            ),
        )
        return response.text


def _make_backend(model_type: str, model_name: str, api_key: str | None) -> _LlmBackend:
    model_type = model_type.lower()
    if model_type == "ollama":
        try:
            import ollama  # noqa: F401
        except ImportError:
            raise ImportError("Package 'ollama' is required. Install: pip install ollama")
        return OllamaBackend(model_name)
    if model_type == "openai":
        if not api_key:
            raise ValueError("api_key required for OpenAI")
        try:
            import openai  # noqa: F401
        except ImportError:
            raise ImportError("Package 'openai' is required. Install: pip install openai")
        return OpenAIBackend(model_name, api_key)
    if model_type == "anthropic":
        if not api_key:
            raise ValueError("api_key required for Anthropic")
        try:
            import anthropic  # noqa: F401
        except ImportError:
            raise ImportError("Package 'anthropic' is required. Install: pip install anthropic")
        return AnthropicBackend(model_name, api_key)
    if model_type == "google":
        if not api_key:
            raise ValueError("api_key required for Google")
        try:
            from google import genai  # noqa: F401
        except ImportError:
            raise ImportError("Package 'google-genai' is required. Install: pip install google-genai")
        return GoogleBackend(model_name, api_key)
    raise ValueError(f"Unknown model_type '{model_type}'. Choose: ollama, openai, anthropic, google")


# ---------------------------------------------------------------------------
# ResponseParser
# ---------------------------------------------------------------------------

_VALID_CLASSIFICATIONS = {"viral-known", "viral-unknown", "non-viral"}

# HTTP / gRPC status codes that indicate a transient server overload and are
# worth retrying (429 = rate-limit, 5xx = server-side errors).
_RETRYABLE_RE = re.compile(r"\b(429|5\d{2})\b")


class ResponseParseError(ValueError):
    """The model replied, but the reply is not a usable assessment."""


class ResponseParser:
    """
    Parses the raw LLM text response into an LlmOutput dataclass.

    The JSON object is located leniently (code fences, <think> blocks and prose
    around it are tolerated — common with local models), but its content is
    validated strictly: a reply without a numeric vq_score or with an unknown
    classification is a parse-error, never silently turned into "non-viral".
    """

    _THINK_RE = re.compile(r"<think>.*?</think>", re.IGNORECASE | re.DOTALL)
    _FENCE_RE = re.compile(r"```(?:json)?", re.IGNORECASE)

    @classmethod
    def _extract_json(cls, raw: str) -> dict:
        text = cls._THINK_RE.sub("", raw or "")
        text = cls._FENCE_RE.sub("", text).strip()
        try:
            data = json.loads(text)
        except json.JSONDecodeError:
            # Fallback: first decodable JSON object inside surrounding prose.
            decoder, data = json.JSONDecoder(), None
            for i, ch in enumerate(text):
                if ch != "{":
                    continue
                try:
                    data, _ = decoder.raw_decode(text[i:])
                    break
                except json.JSONDecodeError:
                    continue
            if data is None:
                raise ResponseParseError("no JSON object in the reply")
        if not isinstance(data, dict):
            raise ResponseParseError("reply JSON is not an object")
        return data

    @staticmethod
    def _normalise_class(value) -> str:
        """'Viral Known' / 'viral_known' → 'viral-known'."""
        return re.sub(r"[\s_]+", "-", str(value or "").strip().lower())

    @classmethod
    def parse(cls, raw: str, seq_id: str, model: str, mode: str) -> LlmOutput:
        data = cls._extract_json(raw)       # raises ResponseParseError

        score = data.get("vq_score")
        if isinstance(score, bool) or score is None:
            raise ResponseParseError("missing vq_score")
        try:
            vq_score = max(0, min(100, int(round(float(score)))))
        except (TypeError, ValueError):
            raise ResponseParseError(f"vq_score is not a number: {score!r}") from None

        classification = cls._normalise_class(data.get("classification"))
        if classification not in _VALID_CLASSIFICATIONS:
            raise ResponseParseError(
                f"unknown classification {data.get('classification')!r}")

        # Novelty is optional: a missing or unknown tier is left empty rather than
        # failing an otherwise valid assessment.
        novelty = cls._normalise_class(data.get("novelty"))
        if novelty not in NOVELTY_TIERS:
            if novelty:
                logger.warning(f"{seq_id}: unknown novelty {data.get('novelty')!r} — ignored.")
            novelty = ""

        return LlmOutput(
            seq_id=seq_id,
            model=model,
            mode=mode,
            vq_score=vq_score,
            classification=classification,
            analysis=str(data.get("analysis", "")).strip(),
            blastn_species=str(data.get("blastn_species", "")).strip(),
            novelty=novelty,
        )

    @staticmethod
    def error_output(seq_id: str, model: str, mode: str, kind: str, message: str) -> LlmOutput:
        """An LlmOutput that records a failure — kind is 'api-error' or 'parse-error'."""
        assert kind in LLM_ERROR_CLASSES
        return LlmOutput(seq_id=seq_id, model=model, mode=mode, vq_score=0,
                         classification=kind, analysis="", error=message)


# ---------------------------------------------------------------------------
# SequenceScorer  (main entry point)
# ---------------------------------------------------------------------------

class SequenceScorer:
    """
    Scores NucSequence objects using an LLM and attaches an LlmOutput to each.

    Parameters
    ----------
    model_type : str
        One of "ollama", "openai", "anthropic", "google".
    model_name : str
        Model identifier (e.g. "llama3.2", "gpt-4o", "claude-sonnet-4-6").
    mode : LlmMode
        HIGH — all hits, full family description.
        LOW  — best hits only, standardised family description.
    api_key : str | None
        Required for openai / anthropic / google.
    low_mode_delay : float
        Seconds to wait between requests when running in LOW mode. Default 30.
        Set to 0 to disable. HIGH mode never waits.
    """

    def __init__(
        self,
        model_type: str,
        model_name: str,
        mode: LlmMode = LlmMode.LOW,
        api_key: str | None = None,
        low_mode_delay: float = 30.0,
    ):
        self._mode            = mode
        self._model_name      = model_name
        self._low_mode_delay  = low_mode_delay
        self._backend         = _make_backend(model_type, model_name, api_key)
        self._system_prompt   = PromptBuilder.system_prompt()

    def score(self, nuc_seqs: list[NucSequence]) -> list[LlmOutput]:
        """
        Scores each sequence, attaches the result to seq.llm_output, and
        returns the full list of LlmOutput objects in the same order.
        The backend is released at the end — even on error or Ctrl+C — so a
        local Ollama model does not linger in memory.
        """
        outputs: list[LlmOutput] = []
        try:
            for i, seq in enumerate(nuc_seqs):
                logger.info(
                    f"[{i+1}/{len(nuc_seqs)}] Scoring {seq.id} "
                    f"({self._mode.value}-token, model={self._model_name})"
                )
                result = self._score_one(seq)
                seq.llm_output = result
                outputs.append(result)

                if result.error:
                    logger.error(f"  ✗ {seq.id} [{result.classification}]: {result.error}")
                else:
                    logger.success(
                        f"  ✓ {seq.id}: score={result.vq_score}, "
                        f"class={result.classification}"
                    )

                if (self._mode == LlmMode.LOW and self._low_mode_delay > 0
                        and self._backend.needs_rate_limit
                        and i < len(nuc_seqs) - 1 and not result.error):
                    logger.debug(f"Low-token mode: waiting {self._low_mode_delay}s before next request...")
                    time.sleep(self._low_mode_delay)
        finally:
            try:
                self._backend.release()
            except Exception as exc:          # never let cleanup mask the real outcome
                logger.warning(f"Could not release the LLM backend: {exc}")

        return outputs

    # Retry constants for transient backend failures.
    _MAX_RETRIES  = 3
    _RETRY_DELAYS = (15.0, 30.0, 60.0)   # seconds between successive attempts
    # A reply that cannot be interpreted is asked for again this many times
    # (models — local ones especially — are non-deterministic at temperature 0.3).
    _PARSE_RETRIES = 1

    def _score_one(self, seq: NucSequence) -> LlmOutput:
        user_json = PromptBuilder.build_input(seq, self._mode)
        last_exc: Exception | None = None
        api_attempt = parse_attempt = 0
        while True:
            try:
                raw = self._backend.call(self._system_prompt, user_json)
            except Exception as exc:
                last_exc = exc
                if _RETRYABLE_RE.search(str(exc)) and api_attempt < self._MAX_RETRIES:
                    delay = self._RETRY_DELAYS[api_attempt]
                    api_attempt += 1
                    logger.warning(
                        f"Backend transient error for {seq.id} "
                        f"(attempt {api_attempt}/{self._MAX_RETRIES}): {exc} — "
                        f"retrying in {delay:.0f}s ..."
                    )
                    time.sleep(delay)
                    continue
                logger.error(f"Backend call failed for {seq.id}: {last_exc}")
                return ResponseParser.error_output(
                    seq.id, self._model_name, self._mode.value, "api-error", str(last_exc))

            try:
                return ResponseParser.parse(raw, seq.id, self._model_name, self._mode.value)
            except ResponseParseError as exc:
                if parse_attempt < self._PARSE_RETRIES:
                    parse_attempt += 1
                    logger.warning(
                        f"Unreadable reply for {seq.id} ({exc}) — asking again "
                        f"({parse_attempt}/{self._PARSE_RETRIES}) ..."
                    )
                    continue
                snippet = " ".join((raw or "").split())[:300]
                return ResponseParser.error_output(
                    seq.id, self._model_name, self._mode.value, "parse-error",
                    f"Invalid model response: {exc}. Reply started with: {snippet!r}")
