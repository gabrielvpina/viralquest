import json
import re
import time
from enum import Enum
from pathlib import Path

from loguru import logger

from viralquest.biodata import (
    BlastnResult,
    BlastxResult,
    HmmDomain,
    LlmOutput,
    NucSequence,
    Orf,
)


class LlmMode(Enum):
    HIGH = "high"   # all hits, Pfam details included, full family description
    LOW  = "low"    # best hits only, Pfam details omitted, standardised description


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
        "subject":       h.stitle,
        "pct_identity":  h.pident,
        "query_coverage": h.qcovhsp,
        "e_value":       h.evalue,
    }


def _domain_dict(d: HmmDomain, include_details: bool) -> dict:
    r: dict = {
        "database":    d.database,
        "target":      d.target,
        "score":       d.score,
        "e_value":     d.e_value,
        "description": d.description,
        "type":        d.type,
    }
    if include_details:
        r["details"] = d.details
    return r


def _orf_dict(orf: Orf, include_details: bool) -> dict:
    return {
        "name":        orf.name,
        "strand":      orf.strand,
        "frame":       orf.frame,
        "length_aa":   orf.length_aa,
        "hmm_domains": [_domain_dict(d, include_details) for d in orf.domains],
    }


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
        include_details = high   # Pfam details only in high-token mode

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

        return {
            "sequence_id":      seq.id,
            "sequence_stats": {
                "length_nt":   seq.length,
                "gc_content":  seq.gc_content,
                "n_count":     seq.n_count,
            },
            "is_viral_flag":    seq.is_viral,
            "taxonomy":         taxonomy_dict,
            "blastx_hits":      blastx_hits,
            "blastn_hits":      blastn_hits,
            "orfs":             [_orf_dict(o, include_details) for o in seq.orfs],
            "viral_family_info": family_info,
        }


# ---------------------------------------------------------------------------
# LLM backends  (one class per provider)
# ---------------------------------------------------------------------------

class _LlmBackend:
    def call(self, system_prompt: str, user_json: dict) -> str:
        raise NotImplementedError


class OllamaBackend(_LlmBackend):
    """Direct Ollama backend. Uses think=False to suppress chain-of-thought."""

    def __init__(self, model_name: str):
        self.model_name = model_name

    def call(self, system_prompt: str, user_json: dict) -> str:
        import ollama
        messages = [
            {"role": "system", "content": system_prompt},
            {"role": "user",   "content": json.dumps(user_json, ensure_ascii=False)},
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
                {"role": "user",   "content": json.dumps(user_json, ensure_ascii=False)},
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
            messages=[{"role": "user", "content": json.dumps(user_json, ensure_ascii=False)}],
        )
        # take only text blocks — skips any extended-thinking blocks automatically
        return "".join(block.text for block in response.content if block.type == "text")


class GoogleBackend(_LlmBackend):
    """Google Gemini backend. Uses json MIME type to enforce structured output."""

    def __init__(self, model_name: str, api_key: str):
        self.model_name = model_name
        self.api_key    = api_key

    def call(self, system_prompt: str, user_json: dict) -> str:
        import google.generativeai as genai
        genai.configure(api_key=self.api_key)
        model = genai.GenerativeModel(
            model_name=self.model_name,
            system_instruction=system_prompt,
            generation_config=genai.GenerationConfig(
                temperature=0.3,
                response_mime_type="application/json",
            ),
        )
        response = model.generate_content(json.dumps(user_json, ensure_ascii=False))
        return response.text


def _make_backend(model_type: str, model_name: str, api_key: str | None) -> _LlmBackend:
    model_type = model_type.lower()
    if model_type == "ollama":
        return OllamaBackend(model_name)
    if model_type == "openai":
        if not api_key:
            raise ValueError("api_key required for OpenAI")
        return OpenAIBackend(model_name, api_key)
    if model_type == "anthropic":
        if not api_key:
            raise ValueError("api_key required for Anthropic")
        return AnthropicBackend(model_name, api_key)
    if model_type == "google":
        if not api_key:
            raise ValueError("api_key required for Google")
        return GoogleBackend(model_name, api_key)
    raise ValueError(f"Unknown model_type '{model_type}'. Choose: ollama, openai, anthropic, google")


# ---------------------------------------------------------------------------
# ResponseParser
# ---------------------------------------------------------------------------

_VALID_CLASSIFICATIONS = {"viral-known", "viral-unknown", "non-viral"}


class ResponseParser:
    """Parses the raw LLM text response into an LlmOutput dataclass."""

    @staticmethod
    def parse(raw: str, seq_id: str, model: str, mode: str) -> LlmOutput:
        try:
            text = raw.strip()
            # strip markdown code fences if the model added them
            text = re.sub(r"^```(?:json)?\s*", "", text, flags=re.IGNORECASE)
            text = re.sub(r"\s*```\s*$",        "", text)
            data = json.loads(text)

            vq_score = max(0, min(100, int(data.get("vq_score", 0))))

            classification = str(data.get("classification", "")).strip().lower()
            if classification not in _VALID_CLASSIFICATIONS:
                logger.warning(
                    f"Unexpected classification '{classification}' for {seq_id}; defaulting to non-viral"
                )
                classification = "non-viral"

            analysis = str(data.get("analysis", "")).strip()

            return LlmOutput(
                seq_id=seq_id,
                model=model,
                mode=mode,
                vq_score=vq_score,
                classification=classification,
                analysis=analysis,
            )

        except Exception as exc:
            return LlmOutput(
                seq_id=seq_id,
                model=model,
                mode=mode,
                vq_score=0,
                classification="non-viral",
                analysis="",
                error=f"ParseError: {exc} | raw[:300]: {raw[:300]}",
            )


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
        HIGH — all hits, Pfam details, full family description.
        LOW  — best hits only, no Pfam details, standardised description.
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
        """
        outputs: list[LlmOutput] = []
        for i, seq in enumerate(nuc_seqs):
            logger.info(
                f"[{i+1}/{len(nuc_seqs)}] Scoring {seq.id} "
                f"({self._mode.value}-token, model={self._model_name})"
            )
            result = self._score_one(seq)
            seq.llm_output = result
            outputs.append(result)

            if result.error:
                logger.error(f"  ✗ {seq.id}: {result.error}")
            else:
                logger.success(
                    f"  ✓ {seq.id}: score={result.vq_score}, "
                    f"class={result.classification}"
                )

            if self._mode == LlmMode.LOW and self._low_mode_delay > 0 and i < len(nuc_seqs) - 1:
                logger.debug(f"Low-token mode: waiting {self._low_mode_delay}s before next request...")
                time.sleep(self._low_mode_delay)

        return outputs

    def _score_one(self, seq: NucSequence) -> LlmOutput:
        try:
            user_json = PromptBuilder.build_input(seq, self._mode)
            raw       = self._backend.call(self._system_prompt, user_json)
            return ResponseParser.parse(raw, seq.id, self._model_name, self._mode.value)
        except Exception as exc:
            logger.error(f"Backend call failed for {seq.id}: {exc}")
            return LlmOutput(
                seq_id=seq.id,
                model=self._model_name,
                mode=self._mode.value,
                vq_score=0,
                classification="non-viral",
                analysis="",
                error=str(exc),
            )
