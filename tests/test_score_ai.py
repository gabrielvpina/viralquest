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
    ResponseParseError,
    ResponseParser,
    SequenceScorer,
    _blastx_dict,
    _blastn_dict,
    _domain_dict,
    _make_backend,
    _orf_dict,
    _to_user_message,
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
        bit_score=200.0,
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
    # Pfam `details` text is never sent to the LLM (any mode).
    def test_details_never_included(self):
        d = _domain_dict(make_domain(details="Pfam details here"))
        assert "details" not in d

    def test_core_fields_always_present(self):
        d = _domain_dict(make_domain())
        assert set(d) == {"database", "target", "score", "e_value", "description", "type"}

    def test_empty_type_omitted(self):
        dom = make_domain(); dom.type = ""          # RVDB / Vfam / EggNOG
        assert "type" not in _domain_dict(dom)

    def test_score_and_evalue_trimmed(self):
        dom = make_domain(score=152.4567); dom.e_value = 1.23456e-45
        d = _domain_dict(dom)
        assert d["score"] == 152.5 and d["e_value"] == 1.2e-45


class TestOrfDict:
    def test_no_aa_or_nuc_sequence(self):
        d = _orf_dict(make_orf())
        assert "aa_sequence" not in d and "nuc_sequence" not in d

    def test_domains_serialised_without_details(self):
        orf = make_orf(domains=[make_domain(details="text")])
        d = _orf_dict(orf)
        assert len(d["hmm_domains"]) == 1
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

    def test_output_schema_fields(self):
        data = json.loads(PromptBuilder.system_prompt())
        schema = data["output"]["schema"]
        assert set(schema) == {"vq_score", "classification", "blastn_species", "novelty", "analysis"}

    def test_novelty_tiers_defined_in_prompt(self):
        from viralquest.biodata import NOVELTY_TIERS
        tiers = json.loads(PromptBuilder.system_prompt())["novelty"]
        assert set(NOVELTY_TIERS) <= set(tiers)

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
        common = {"sequence_id", "sequence_stats", "is_viral_flag", "taxonomy",
                  "blastx_hits", "blastn_hits", "unannotated_orfs", "viral_family_info"}
        assert (common | {"orfs"}).issubset(PromptBuilder.build_input(make_seq(), LlmMode.HIGH))
        assert (common | {"orf_summary"}).issubset(PromptBuilder.build_input(make_seq(), LlmMode.LOW))

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

    def test_high_mode_pfam_details_omitted(self):
        seq = make_seq()
        seq.orfs.append(make_orf(domains=[make_domain(database="Pfam", details="Pfam detail text")]))
        d = PromptBuilder.build_input(seq, LlmMode.HIGH)
        assert "details" not in d["orfs"][0]["hmm_domains"][0]

    def test_low_mode_pfam_details_omitted(self):
        seq = make_seq()
        seq.orfs.append(make_orf(domains=[make_domain(database="Pfam", details="Pfam detail text")]))
        d = PromptBuilder.build_input(seq, LlmMode.LOW)
        assert "Pfam detail text" not in json.dumps(d)

    def test_low_mode_non_pfam_details_also_omitted(self):
        seq = make_seq()
        seq.orfs.append(make_orf(domains=[make_domain(database="RVDB", details="some detail")]))
        d = PromptBuilder.build_input(seq, LlmMode.LOW)
        assert "some detail" not in json.dumps(d)

    # --- ORF payload: annotated vs unannotated ---

    def _seq_with_orfs(self):
        seq = make_seq()
        seq.orfs += [
            make_orf("seq1_ORF_1_900+", domains=[make_domain(database="Pfam", score=152.4)]),
            make_orf("seq1_ORF_950_1100+"),
            make_orf("seq1_ORF_1200_1400-"),
            make_orf("seq1_ORF_1500_1600+"),
            make_orf("seq1_ORF_1700_1800+"),
        ]
        for o, aa in zip(seq.orfs[1:], (412, 60, 260, 198)):
            o.length_aa = aa
        return seq

    def test_high_lists_only_annotated_orfs(self):
        d = PromptBuilder.build_input(self._seq_with_orfs(), LlmMode.HIGH)
        assert [o["name"] for o in d["orfs"]] == ["seq1_ORF_1_900+"]
        assert d["orfs"][0]["hmm_domains"][0]["target"] == "PF00001"

    def test_unannotated_orfs_summarised_in_both_modes(self):
        for mode in (LlmMode.HIGH, LlmMode.LOW):
            d = PromptBuilder.build_input(self._seq_with_orfs(), mode)
            assert d["unannotated_orfs"] == {"count": 4, "longest_aa": [412, 260, 198]}

    def test_low_summarises_annotated_orfs_with_domain_names(self):
        d = PromptBuilder.build_input(self._seq_with_orfs(), LlmMode.LOW)
        summ = d["orf_summary"]
        assert summ["total"] == 5 and summ["with_domains"] == 1
        # seq-id prefix dropped; domain name, description and rounded score kept
        assert summ["annotated_orfs"] == [
            "ORF_1_900+ (100 aa): Pfam PF00001 – RNA-dep RNA polymerase [152]"]
        assert "orfs" not in d

    def test_low_lists_every_domain_of_an_orf(self):
        seq = make_seq()
        seq.orfs.append(make_orf("seq1_ORF_1_900+", domains=[
            make_domain(database="RVDB", score=98.0), make_domain(database="Pfam", score=150.0)]))
        line = PromptBuilder.build_input(seq, LlmMode.LOW)["orf_summary"]["annotated_orfs"][0]
        assert "RVDB PF00001" in line and "Pfam PF00001" in line and line.count(";") == 1

    def test_no_orfs_gives_empty_summaries(self):
        d = PromptBuilder.build_input(make_seq(), LlmMode.LOW)
        assert d["orf_summary"] == {"total": 0, "with_domains": 0, "annotated_orfs": []}
        assert d["unannotated_orfs"] == {"count": 0, "longest_aa": []}


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

    def test_unknown_classification_is_a_parse_error(self):
        # Never silently turned into "non-viral".
        raw = json.dumps({"vq_score": 50, "classification": "maybe-viral", "analysis": "x"})
        with pytest.raises(ResponseParseError, match="classification"):
            self._parse(raw)

    def test_classification_spaces_and_underscores_normalised(self):
        for variant in ("Viral Known", "viral_known", " VIRAL-KNOWN "):
            raw = json.dumps({"vq_score": 90, "classification": variant})
            assert self._parse(raw).classification == "viral-known"

    def test_classification_normalised_to_lowercase(self):
        raw = json.dumps({"vq_score": 50, "classification": "Viral-Known", "analysis": "x"})
        r = self._parse(raw)
        assert r.classification == "viral-known"

    def test_invalid_json_is_a_parse_error(self):
        with pytest.raises(ResponseParseError):
            self._parse("this is not json")

    def test_empty_string_is_a_parse_error(self):
        with pytest.raises(ResponseParseError):
            self._parse("")

    def test_missing_vq_score_is_a_parse_error(self):
        raw = json.dumps({"classification": "non-viral", "analysis": "x"})
        with pytest.raises(ResponseParseError, match="vq_score"):
            self._parse(raw)

    def test_non_numeric_vq_score_is_a_parse_error(self):
        raw = json.dumps({"vq_score": "high", "classification": "non-viral"})
        with pytest.raises(ResponseParseError, match="vq_score"):
            self._parse(raw)

    def test_numeric_string_score_accepted(self):
        raw = json.dumps({"vq_score": "72.6", "classification": "viral-unknown"})
        assert self._parse(raw).vq_score == 73

    def test_json_array_is_a_parse_error(self):
        with pytest.raises(ResponseParseError):
            self._parse("[1, 2, 3]")

    # --- lenient extraction (common with local models) ---

    def test_json_surrounded_by_prose_is_extracted(self):
        raw = "Sure! Here is my assessment:\n" + good_response(64, "viral-unknown") + "\nHope it helps."
        r = self._parse(raw)
        assert r.vq_score == 64 and r.classification == "viral-unknown"

    def test_think_block_is_ignored(self):
        raw = '<think>maybe {"vq_score": 1} hmm</think>\n' + good_response(88, "viral-known")
        assert self._parse(raw).vq_score == 88

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

    def test_novelty_parsed(self):
        raw = json.dumps({"vq_score": 70, "classification": "viral-unknown", "novelty": "Novel Species"})
        assert self._parse(raw).novelty == "novel-species"

    def test_novelty_missing_or_unknown_is_empty_not_an_error(self):
        assert self._parse(json.dumps({"vq_score": 70, "classification": "viral-unknown"})).novelty == ""
        raw = json.dumps({"vq_score": 70, "classification": "viral-unknown", "novelty": "weird"})
        out = self._parse(raw)
        assert out.novelty == "" and out.classification == "viral-unknown"

    def test_error_output_marks_the_failure(self):
        r = ResponseParser.error_output("s1", "m", "low", "parse-error", "bad reply")
        assert r.classification == "parse-error" and r.error == "bad reply"
        assert r.vq_score == 0 and r.blastn_species == "" and r.analysis == ""

    def test_error_output_rejects_unknown_kind(self):
        with pytest.raises(AssertionError):
            ResponseParser.error_output("s1", "m", "low", "non-viral", "x")

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


# ---------------------------------------------------------------------------
# Compact user message
# ---------------------------------------------------------------------------

class TestUserMessage:
    def test_compact_separators(self):
        msg = _to_user_message({"a": 1, "b": [1, 2]})
        assert msg == '{"a":1,"b":[1,2]}'

    def test_roundtrips_and_keeps_unicode(self):
        payload = {"desc": "RNA-dependent – polymerase", "n": 3}
        msg = _to_user_message(payload)
        assert json.loads(msg) == payload and "–" in msg

    def test_every_backend_sends_the_compact_message(self):
        import inspect
        from viralquest import score_ai
        src = inspect.getsource(score_ai)
        assert src.count("_to_user_message(user_json)") == 4
        assert src.count("json.dumps(user_json") == 1   # only inside the helper itself


# ---------------------------------------------------------------------------
# SequenceScorer — invalid replies and backend release
# ---------------------------------------------------------------------------

class TestScorerInvalidReplies:
    def _scorer(self, replies):
        s = SequenceScorer("ollama", "llama3", mode=LlmMode.LOW, low_mode_delay=0.0)
        s._backend = MagicMock()
        s._backend.call.side_effect = replies
        return s

    def test_invalid_reply_retried_once_then_succeeds(self):
        s = self._scorer(["not json", good_response(70, "viral-unknown")])
        out = s.score([make_seq()])[0]
        assert out.classification == "viral-unknown" and out.error is None
        assert s._backend.call.call_count == 2

    def test_invalid_twice_becomes_parse_error_not_non_viral(self):
        s = self._scorer(["not json", "still not json"])
        out = s.score([make_seq()])[0]
        assert out.classification == "parse-error"
        assert "Invalid model response" in out.error and "still not json" in out.error
        assert s._backend.call.call_count == 2

    def test_backend_exception_is_api_error(self):
        s = self._scorer([RuntimeError("connection refused")])
        out = s.score([make_seq()])[0]
        assert out.classification == "api-error" and "connection refused" in out.error

    def test_one_bad_sequence_does_not_stop_the_others(self):
        s = self._scorer(["x", "y", good_response(91, "viral-known")])
        outs = s.score([make_seq("a"), make_seq("b")])
        assert [o.classification for o in outs] == ["parse-error", "viral-known"]


class TestBackendRelease:
    def _scorer(self):
        s = SequenceScorer("ollama", "llama3", mode=LlmMode.HIGH)
        s._backend = MagicMock()
        return s

    def test_released_after_scoring(self):
        s = self._scorer()
        s._backend.call.return_value = good_response()
        s.score([make_seq()])
        s._backend.release.assert_called_once()

    def test_released_even_when_scoring_is_interrupted(self):
        s = self._scorer()
        s._backend.call.side_effect = KeyboardInterrupt
        with pytest.raises(KeyboardInterrupt):
            s.score([make_seq()])
        s._backend.release.assert_called_once()

    def test_release_failure_does_not_break_scoring(self):
        s = self._scorer()
        s._backend.call.return_value = good_response(80, "viral-known")
        s._backend.release.side_effect = ConnectionError("ollama down")
        assert s.score([make_seq()])[0].vq_score == 80

    def test_ollama_release_unloads_with_keep_alive_zero(self):
        import ollama
        with patch.object(ollama, "generate") as gen:
            OllamaBackend("qwen3:4b").release()
        gen.assert_called_once_with(model="qwen3:4b", keep_alive=0)

    def test_cloud_backends_release_is_a_no_op(self):
        from viralquest.score_ai import _LlmBackend
        assert _LlmBackend().release() is None
