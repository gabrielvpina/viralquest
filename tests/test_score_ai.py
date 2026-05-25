import json
import time
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from viralquest.biodata import (
    BlastnResult,
    BlastxResult,
    HmmDomain,
    LlmOutput,
    NucSequence,
    Orf,
    Taxonomy,
    ViralFamilyInfo,
)
from viralquest.score_ai import (
    LlmMode,
    PromptBuilder,
    ResponseParser,
    SequenceScorer,
    _blastx_dict,
    _blastn_dict,
    _domain_dict,
    _make_backend,
    _orf_dict,
    AnthropicBackend,
    GoogleBackend,
    OllamaBackend,
    OpenAIBackend,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_seq(seq_id: str = "seq1") -> NucSequence:
    return NucSequence(id=seq_id, sequence="ATGCGTNNNATG")


def make_blastx(subject_title: str = "RdRp [Influenza A virus]", bit_score: float = 300.0) -> BlastxResult:
    h = BlastxResult(
        query_id="seq1", subject_id="NC_001",
        subject_title=subject_title,
        pct_identity=95.0, aln_length=200, mismatches=5,
        gap_opens=0, query_start=1, query_end=200,
        subject_start=1, subject_end=200,
        e_value=1e-80, bit_score=bit_score,
    )
    h.query_coverage = 85.0
    return h


def make_blastn(stitle: str = "Influenza A virus genome") -> BlastnResult:
    return BlastnResult(
        qseqid="seq1", qlen=1200, slen=13600,
        qcovhsp=80, pident=92.5, evalue=1e-60,
        stitle=stitle,
    )


def make_domain(database: str = "RVDB", score: float = 120.0, details: str = "long Pfam detail text") -> HmmDomain:
    return HmmDomain(
        database=database, target="PF00001",
        score=score, e_value=1e-10,
        start=1, stop=50, length=49,
        description="RNA-dep RNA polymerase",
        type="Family",
        details=details,
    )


def make_orf(name: str = "orf1", domains: list | None = None) -> Orf:
    orf = Orf(
        start_codon="ATG", stop_codon="TAA",
        start_position=1, stop_position=300,
        strand="+", frame=1,
        length_aa=100, length_nt=300,
        bigger_than_50=True,
        aa_sequence="MKVLS",
        nuc_sequence="ATGAAAGTTCTGTCT",
        orf_type="complete",
        name=name,
    )
    if domains:
        orf.domains.extend(domains)
    return orf


def make_taxonomy() -> Taxonomy:
    return Taxonomy(
        tax_id=11520,
        scientific_name="Influenza A virus",
        no_rank=None,
        clade="Riboviria",
        kingdom="Orthornavirae",
        phylum="Negarnaviricota",
        class_="Insthoviricetes",
        order="Articulavirales",
        family="Orthomyxoviridae",
        subfamily=None,
        genus="Alphainfluenzavirus",
        species="Alphainfluenzavirus influenzae",
        genome="ssRNA(-)",
    )


def make_vfi(high: str = "Full ICTV text", low: str = "Compact text") -> ViralFamilyInfo:
    return ViralFamilyInfo(
        source="ICTV",
        type="Family",
        name="Orthomyxoviridae",
        info_high=high,
        info_low=low,
    )


def good_response(
    vq: int = 85,
    cls: str = "viral-known",
    analysis: str = "Strong hit.",
    blastn_species: str = "Influenza A virus",
) -> str:
    return json.dumps({
        "vq_score": vq,
        "classification": cls,
        "analysis": analysis,
        "blastn_species": blastn_species,
    })


# ---------------------------------------------------------------------------
# LlmMode
# ---------------------------------------------------------------------------

class TestLlmMode:
    def test_values(self):
        assert LlmMode.HIGH.value == "high"
        assert LlmMode.LOW.value  == "low"


# ---------------------------------------------------------------------------
# Serialisation helpers
# ---------------------------------------------------------------------------

class TestBlastxDict:
    def test_fields_present(self):
        d = _blastx_dict(make_blastx())
        assert set(d) == {"subject", "species", "pct_identity", "query_coverage",
                          "aln_length", "e_value", "bit_score"}

    def test_species_extracted(self):
        d = _blastx_dict(make_blastx("coat protein [Tobacco mosaic virus]"))
        assert d["species"] == "Tobacco mosaic virus"

    def test_no_sequence_field(self):
        d = _blastx_dict(make_blastx())
        assert "sequence" not in d and "aa_sequence" not in d


class TestBlastnDict:
    def test_fields_present(self):
        d = _blastn_dict(make_blastn())
        assert set(d) == {"subject", "pct_identity", "query_coverage", "e_value"}

    def test_no_sequence_field(self):
        assert "sequence" not in _blastn_dict(make_blastn())


class TestDomainDict:
    def test_details_included_when_flag_true(self):
        d = _domain_dict(make_domain(details="Pfam details here"), include_details=True)
        assert d["details"] == "Pfam details here"

    def test_details_omitted_when_flag_false(self):
        d = _domain_dict(make_domain(details="Pfam details here"), include_details=False)
        assert "details" not in d

    def test_core_fields_always_present(self):
        d = _domain_dict(make_domain(), include_details=False)
        assert {"database", "target", "score", "e_value", "description", "type"}.issubset(d)


class TestOrfDict:
    def test_no_aa_or_nuc_sequence(self):
        d = _orf_dict(make_orf(), include_details=False)
        assert "aa_sequence" not in d and "nuc_sequence" not in d

    def test_domains_serialised(self):
        orf = make_orf(domains=[make_domain()])
        d = _orf_dict(orf, include_details=True)
        assert len(d["hmm_domains"]) == 1
        assert "details" in d["hmm_domains"][0]

    def test_domain_details_omitted_in_low(self):
        orf = make_orf(domains=[make_domain()])
        d = _orf_dict(orf, include_details=False)
        assert "details" not in d["hmm_domains"][0]


# ---------------------------------------------------------------------------
# PromptBuilder.system_prompt
# ---------------------------------------------------------------------------

class TestPromptBuilderSystemPrompt:
    def test_loads_valid_json(self):
        sp = PromptBuilder.system_prompt()
        data = json.loads(sp)
        assert isinstance(data, dict)

    def test_required_keys_present(self):
        data = json.loads(PromptBuilder.system_prompt())
        assert {"role", "task", "scoring", "classification", "output"}.issubset(data)

    def test_output_schema_has_four_fields(self):
        data = json.loads(PromptBuilder.system_prompt())
        schema = data["output"]["schema"]
        assert set(schema) == {"vq_score", "classification", "blastn_species", "analysis"}

    def test_field_extraction_blastn_species_defined(self):
        data = json.loads(PromptBuilder.system_prompt())
        assert "field_extraction" in data
        assert "blastn_species" in data["field_extraction"]

    def test_all_three_classifications_defined(self):
        data = json.loads(PromptBuilder.system_prompt())
        assert set(data["classification"]) == {"viral-known", "viral-unknown", "non-viral"}

    def test_same_prompt_for_both_modes(self):
        # system prompt is mode-independent; data volume is controlled in build_input
        sp = PromptBuilder.system_prompt()
        assert isinstance(sp, str) and len(sp) > 0


# ---------------------------------------------------------------------------
# PromptBuilder.build_input
# ---------------------------------------------------------------------------

class TestPromptBuilderBuildInput:
    def test_required_top_level_keys(self):
        seq = make_seq()
        d = PromptBuilder.build_input(seq, LlmMode.LOW)
        assert {"sequence_id", "sequence_stats", "is_viral_flag",
                "taxonomy", "blastx_hits", "blastn_hits",
                "orfs", "viral_family_info"}.issubset(d)

    def test_sequences_never_in_output(self):
        seq = make_seq()
        seq.orfs.append(make_orf())
        for mode in (LlmMode.HIGH, LlmMode.LOW):
            d = PromptBuilder.build_input(seq, mode)
            text = json.dumps(d)
            assert "aa_sequence"  not in text
            assert "nuc_sequence" not in text

    def test_sequence_stats_values(self):
        seq = make_seq()
        d = PromptBuilder.build_input(seq, LlmMode.LOW)
        stats = d["sequence_stats"]
        assert stats["length_nt"] == seq.length
        assert stats["gc_content"] == seq.gc_content

    # --- taxonomy ---

    def test_no_taxonomy_gives_none(self):
        seq = make_seq()
        assert PromptBuilder.build_input(seq, LlmMode.LOW)["taxonomy"] is None

    def test_high_mode_taxonomy_includes_full_lineage(self):
        seq = make_seq()
        seq.taxonomy = make_taxonomy()
        d = PromptBuilder.build_input(seq, LlmMode.HIGH)["taxonomy"]
        assert {"scientific_name", "clade", "kingdom", "phylum",
                "class", "order", "family", "genus", "species", "genome"}.issubset(d)

    def test_low_mode_taxonomy_is_compact(self):
        seq = make_seq()
        seq.taxonomy = make_taxonomy()
        d = PromptBuilder.build_input(seq, LlmMode.LOW)["taxonomy"]
        assert set(d) == {"family", "genus", "order", "genome"}
        assert "scientific_name" not in d

    # --- viral family info ---

    def test_no_vfi_gives_none(self):
        seq = make_seq()
        assert PromptBuilder.build_input(seq, LlmMode.HIGH)["viral_family_info"] is None

    def test_high_mode_uses_info_high(self):
        seq = make_seq()
        seq.viral_family_info = make_vfi(high="FULL TEXT", low="compact")
        d = PromptBuilder.build_input(seq, LlmMode.HIGH)["viral_family_info"]
        assert d["info"] == "FULL TEXT"

    def test_low_mode_uses_info_low(self):
        seq = make_seq()
        seq.viral_family_info = make_vfi(high="FULL TEXT", low="compact")
        d = PromptBuilder.build_input(seq, LlmMode.LOW)["viral_family_info"]
        assert d["info"] == "compact"

    # --- BLAST hits ---

    def test_high_mode_includes_all_blastx_hits(self):
        seq = make_seq()
        seq.blastx_hits.append(make_blastx(bit_score=100.0))
        seq.blastx_hits.append(make_blastx(bit_score=200.0))
        seq.blastx_hits.append(make_blastx(bit_score=300.0))
        d = PromptBuilder.build_input(seq, LlmMode.HIGH)
        assert len(d["blastx_hits"]) == 3

    def test_low_mode_keeps_only_best_blastx(self):
        seq = make_seq()
        seq.blastx_hits.append(make_blastx(bit_score=100.0))
        seq.blastx_hits.append(make_blastx(bit_score=400.0))
        d = PromptBuilder.build_input(seq, LlmMode.LOW)
        assert len(d["blastx_hits"]) == 1
        assert d["blastx_hits"][0]["bit_score"] == 400.0

    def test_high_mode_includes_all_blastn_hits(self):
        seq = make_seq()
        seq.blastn_hits.append(make_blastn("hit 1"))
        seq.blastn_hits.append(make_blastn("hit 2"))
        d = PromptBuilder.build_input(seq, LlmMode.HIGH)
        assert len(d["blastn_hits"]) == 2

    def test_low_mode_keeps_only_best_blastn(self):
        seq = make_seq()
        seq.blastn_hits.append(make_blastn("hit 1"))
        seq.blastn_hits.append(make_blastn("hit 2"))
        d = PromptBuilder.build_input(seq, LlmMode.LOW)
        assert len(d["blastn_hits"]) == 1

    def test_no_hits_gives_empty_lists(self):
        seq = make_seq()
        d = PromptBuilder.build_input(seq, LlmMode.HIGH)
        assert d["blastx_hits"] == []
        assert d["blastn_hits"] == []

    def test_nr_hits_preferred_over_refseq_hits(self):
        seq = make_seq()
        seq.blastx_hits.append(make_blastx("refseq hit [Virus A]", bit_score=100.0))
        seq.blastx_nr_hits.append(make_blastx("nr hit [Virus B]", bit_score=50.0))
        d = PromptBuilder.build_input(seq, LlmMode.HIGH)
        assert all("Virus B" in h["subject"] for h in d["blastx_hits"])

    # --- HMM domains / Pfam details ---

    def test_high_mode_pfam_details_included(self):
        seq = make_seq()
        seq.orfs.append(make_orf(domains=[make_domain(database="Pfam", details="Pfam detail text")]))
        d = PromptBuilder.build_input(seq, LlmMode.HIGH)
        assert d["orfs"][0]["hmm_domains"][0]["details"] == "Pfam detail text"

    def test_low_mode_pfam_details_omitted(self):
        seq = make_seq()
        seq.orfs.append(make_orf(domains=[make_domain(database="Pfam", details="Pfam detail text")]))
        d = PromptBuilder.build_input(seq, LlmMode.LOW)
        assert "details" not in d["orfs"][0]["hmm_domains"][0]

    def test_low_mode_non_pfam_details_also_omitted(self):
        seq = make_seq()
        seq.orfs.append(make_orf(domains=[make_domain(database="RVDB", details="some detail")]))
        d = PromptBuilder.build_input(seq, LlmMode.LOW)
        assert "details" not in d["orfs"][0]["hmm_domains"][0]


# ---------------------------------------------------------------------------
# ResponseParser
# ---------------------------------------------------------------------------

class TestResponseParser:
    def _parse(self, raw, seq_id="seq1", model="m", mode="low"):
        return ResponseParser.parse(raw, seq_id, model, mode)

    def test_valid_json_parsed_correctly(self):
        r = self._parse(good_response(85, "viral-known", "Good evidence.", "Influenza A virus"))
        assert r.vq_score        == 85
        assert r.classification  == "viral-known"
        assert r.analysis        == "Good evidence."
        assert r.blastn_species  == "Influenza A virus"
        assert r.error is None

    def test_all_three_classifications_accepted(self):
        for cls in ("viral-known", "viral-unknown", "non-viral"):
            r = self._parse(good_response(cls=cls))
            assert r.classification == cls

    def test_score_clamped_below_zero(self):
        r = self._parse(good_response(vq=-10))
        assert r.vq_score == 0

    def test_score_clamped_above_100(self):
        r = self._parse(good_response(vq=150))
        assert r.vq_score == 100

    def test_markdown_fence_stripped(self):
        raw = "```json\n" + good_response(72, "viral-unknown", "Weak.") + "\n```"
        r = self._parse(raw)
        assert r.vq_score == 72
        assert r.error is None

    def test_markdown_fence_without_lang_stripped(self):
        raw = "```\n" + good_response() + "\n```"
        assert self._parse(raw).error is None

    def test_unknown_classification_defaults_to_non_viral(self):
        raw = json.dumps({"vq_score": 50, "classification": "maybe-viral", "analysis": "x"})
        r = self._parse(raw)
        assert r.classification == "non-viral"

    def test_classification_normalised_to_lowercase(self):
        raw = json.dumps({"vq_score": 50, "classification": "Viral-Known", "analysis": "x"})
        r = self._parse(raw)
        assert r.classification == "viral-known"

    def test_invalid_json_returns_error(self):
        r = self._parse("this is not json")
        assert r.error is not None
        assert r.vq_score == 0
        assert r.classification == "non-viral"

    def test_empty_string_returns_error(self):
        r = self._parse("")
        assert r.error is not None

    def test_missing_vq_score_defaults_to_zero(self):
        raw = json.dumps({"classification": "non-viral", "analysis": "x"})
        r = self._parse(raw)
        assert r.vq_score == 0

    def test_blastn_species_extracted(self):
        raw = json.dumps({"vq_score": 85, "classification": "viral-known",
                          "analysis": "x", "blastn_species": "Tobacco mosaic virus"})
        assert self._parse(raw).blastn_species == "Tobacco mosaic virus"

    def test_blastn_species_missing_defaults_to_empty(self):
        raw = json.dumps({"vq_score": 85, "classification": "viral-known", "analysis": "x"})
        assert self._parse(raw).blastn_species == ""

    def test_blastn_species_whitespace_stripped(self):
        raw = json.dumps({"vq_score": 85, "classification": "viral-known",
                          "analysis": "x", "blastn_species": "  Tomato spotted wilt virus  "})
        assert self._parse(raw).blastn_species == "Tomato spotted wilt virus"

    def test_error_path_blastn_species_empty(self):
        r = self._parse("not json at all")
        assert r.blastn_species == ""
        assert r.error is not None

    def test_metadata_fields_set(self):
        r = ResponseParser.parse(good_response(), "myseq", "gpt-4o", "high")
        assert r.seq_id == "myseq"
        assert r.model  == "gpt-4o"
        assert r.mode   == "high"


# ---------------------------------------------------------------------------
# _make_backend
# ---------------------------------------------------------------------------

class TestMakeBackend:
    def test_ollama_no_key_required(self):
        b = _make_backend("ollama", "llama3", None)
        assert isinstance(b, OllamaBackend)

    def test_openai_requires_key(self):
        with pytest.raises(ValueError, match="api_key"):
            _make_backend("openai", "gpt-4o", None)

    def test_openai_created_with_key(self):
        assert isinstance(_make_backend("openai", "gpt-4o", "sk-test"), OpenAIBackend)

    def test_anthropic_requires_key(self):
        with pytest.raises(ValueError, match="api_key"):
            _make_backend("anthropic", "claude-sonnet-4-6", None)

    def test_anthropic_created_with_key(self):
        assert isinstance(_make_backend("anthropic", "claude-sonnet-4-6", "sk-ant"), AnthropicBackend)

    def test_google_requires_key(self):
        with pytest.raises(ValueError, match="api_key"):
            _make_backend("google", "gemini-2.0-flash", None)

    def test_google_created_with_key(self):
        assert isinstance(_make_backend("google", "gemini-2.0-flash", "AIza"), GoogleBackend)

    def test_unknown_type_raises(self):
        with pytest.raises(ValueError, match="Unknown model_type"):
            _make_backend("groq", "llama3", None)

    def test_case_insensitive(self):
        assert isinstance(_make_backend("Ollama", "llama3", None), OllamaBackend)
        assert isinstance(_make_backend("OPENAI", "gpt-4o", "key"), OpenAIBackend)


# ---------------------------------------------------------------------------
# SequenceScorer  (backend mocked — no real LLM calls)
# ---------------------------------------------------------------------------

class TestSequenceScorer:
    @pytest.fixture
    def scorer_low(self):
        s = SequenceScorer("ollama", "llama3", mode=LlmMode.LOW, low_mode_delay=0.0)
        s._backend = MagicMock()
        s._backend.call.return_value = good_response(80, "viral-unknown", "Moderate evidence.")
        return s

    @pytest.fixture
    def scorer_high(self):
        s = SequenceScorer("ollama", "llama3", mode=LlmMode.HIGH, low_mode_delay=0.0)
        s._backend = MagicMock()
        s._backend.call.return_value = good_response(95, "viral-known", "Strong evidence.")
        return s

    def test_score_attaches_llm_output_to_seq(self, scorer_low):
        seq = make_seq()
        scorer_low.score([seq])
        assert seq.llm_output is not None
        assert seq.llm_output.vq_score == 80

    def test_score_returns_list_of_outputs(self, scorer_low):
        seqs = [make_seq("s1"), make_seq("s2")]
        outputs = scorer_low.score(seqs)
        assert len(outputs) == 2
        assert outputs[0].seq_id == "s1"
        assert outputs[1].seq_id == "s2"

    def test_score_calls_backend_once_per_seq(self, scorer_low):
        seqs = [make_seq("s1"), make_seq("s2"), make_seq("s3")]
        scorer_low.score(seqs)
        assert scorer_low._backend.call.call_count == 3

    def test_backend_call_receives_system_prompt_and_json(self, scorer_low):
        seq = make_seq()
        scorer_low.score([seq])
        args = scorer_low._backend.call.call_args
        system_arg, user_arg = args[0]
        assert isinstance(system_arg, str)
        assert json.loads(system_arg)  # system prompt is valid JSON
        assert isinstance(user_arg, dict)
        assert user_arg["sequence_id"] == "seq1"

    def test_backend_exception_stored_as_error(self):
        scorer = SequenceScorer("ollama", "llama3", mode=LlmMode.HIGH, low_mode_delay=0.0)
        scorer._backend = MagicMock()
        scorer._backend.call.side_effect = RuntimeError("connection refused")
        seq = make_seq()
        scorer.score([seq])
        assert seq.llm_output.error is not None
        assert "connection refused" in seq.llm_output.error

    def test_high_mode_output_fields(self, scorer_high):
        seq = make_seq()
        outputs = scorer_high.score([seq])
        assert outputs[0].mode           == "high"
        assert outputs[0].classification == "viral-known"
        assert outputs[0].vq_score       == 95

    def test_low_mode_output_fields(self, scorer_low):
        seq = make_seq()
        outputs = scorer_low.score([seq])
        assert outputs[0].mode == "low"

    def test_empty_list_returns_empty(self, scorer_low):
        assert scorer_low.score([]) == []

    def test_default_low_mode_delay_is_30(self):
        s = SequenceScorer("ollama", "llama3")
        assert s._low_mode_delay == 30.0

    def test_custom_low_mode_delay_accepted(self):
        s = SequenceScorer("ollama", "llama3", low_mode_delay=60.0)
        assert s._low_mode_delay == 60.0

    def test_low_mode_delay_applied_between_requests(self):
        scorer = SequenceScorer("ollama", "llama3", mode=LlmMode.LOW, low_mode_delay=0.5)
        scorer._backend = MagicMock()
        scorer._backend.call.return_value = good_response()
        seqs = [make_seq("s1"), make_seq("s2")]
        t0 = time.monotonic()
        scorer.score(seqs)
        elapsed = time.monotonic() - t0
        assert elapsed >= 0.5

    def test_low_mode_no_delay_after_last_seq(self):
        scorer = SequenceScorer("ollama", "llama3", mode=LlmMode.LOW, low_mode_delay=5.0)
        scorer._backend = MagicMock()
        scorer._backend.call.return_value = good_response()
        t0 = time.monotonic()
        scorer.score([make_seq()])   # single sequence — no delay after last
        assert time.monotonic() - t0 < 1.0

    def test_high_mode_never_waits(self):
        scorer = SequenceScorer("ollama", "llama3", mode=LlmMode.HIGH, low_mode_delay=30.0)
        scorer._backend = MagicMock()
        scorer._backend.call.return_value = good_response()
        seqs = [make_seq("s1"), make_seq("s2")]
        t0 = time.monotonic()
        scorer.score(seqs)
        assert time.monotonic() - t0 < 1.0   # HIGH never sleeps

    def test_model_name_recorded_in_output(self, scorer_low):
        seq = make_seq()
        scorer_low.score([seq])
        assert seq.llm_output.model == "llama3"

    def test_system_prompt_is_valid_json(self):
        s = SequenceScorer("ollama", "llama3")
        data = json.loads(s._system_prompt)
        assert "output" in data
