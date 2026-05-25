import json
from unittest.mock import MagicMock, patch

import pytest

from viralquest.biodata import HmmDomain, NucSequence, Orf
from viralquest.hmm import (
    HMM_ROLES,
    HmmMetadataLoader,
    HmmResultAttacher,
    HmmRole,
    HmmSearcher,
    HmmSequencePreparer,
    HmmViralFlagSetter,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_orf(name: str, aa_sequence: str = "MKVLS") -> Orf:
    return Orf(
        start_codon="ATG",
        stop_codon="TAA",
        start_position=1,
        stop_position=18,
        strand="+",
        frame=1,
        length_aa=len(aa_sequence),
        length_nt=len(aa_sequence) * 3,
        bigger_than_50=False,
        aa_sequence=aa_sequence,
        nuc_sequence="ATGAAAGTTCTGTCT",
        orf_type="complete",
        name=name,
    )


def make_nuc_seq(seq_id: str = "seq1", sequence: str = "ATGCGT") -> NucSequence:
    return NucSequence(id=seq_id, sequence=sequence)


def make_domain(database: str = "RVDB") -> HmmDomain:
    return HmmDomain(
        database=database,
        target="HMM001",
        score=100.0,
        e_value=1e-10,
        start=1,
        stop=50,
        length=49,
        description="",
        type="",
        details="",
    )


# ---------------------------------------------------------------------------
# HmmRole
# ---------------------------------------------------------------------------

class TestHmmRole:
    def test_filter_value(self):
        assert HmmRole.FILTER.value == "filter"

    def test_characterize_value(self):
        assert HmmRole.CHARACTERIZE.value == "characterize"


# ---------------------------------------------------------------------------
# HMM_ROLES
# ---------------------------------------------------------------------------

class TestHmmRoles:
    def test_rvdb_is_filter(self):
        assert HMM_ROLES["RVDB"] == HmmRole.FILTER

    def test_vfam_is_filter(self):
        assert HMM_ROLES["Vfam"] == HmmRole.FILTER

    def test_eggnon_is_filter(self):
        assert HMM_ROLES["EggNOG"] == HmmRole.FILTER

    def test_pfam_is_characterize(self):
        assert HMM_ROLES["Pfam"] == HmmRole.CHARACTERIZE


# ---------------------------------------------------------------------------
# HmmMetadataLoader
# ---------------------------------------------------------------------------

class TestHmmMetadataLoader:
    def test_load_eggnog(self, tmp_path):
        data = [
            {"EggNOG_TargetID": "COG001", "EggNOG_Description": "Viral replicase"},
            {"EggNOG_TargetID": "COG002", "EggNOG_Description": "Capsid protein"},
        ]
        json_file = tmp_path / "eggnog.json"
        json_file.write_text(json.dumps(data))
        result = HmmMetadataLoader.load(str(json_file), "EggNOG")
        assert len(result) == 2
        assert result["COG001"]["description"] == "Viral replicase"
        assert result["COG001"]["type"] == ""
        assert result["COG001"]["details"] == ""

    def test_load_pfam(self, tmp_path):
        data = [
            {
                "Pfam_TargetID": "PF00001",
                "Pfam_Description": "7tm_1",
                "Pfam_Type": "Family",
                "Pfam_Details": "7 transmembrane receptor",
            }
        ]
        json_file = tmp_path / "pfam.json"
        json_file.write_text(json.dumps(data))
        result = HmmMetadataLoader.load(str(json_file), "Pfam")
        assert result["PF00001"]["description"] == "7tm_1"
        assert result["PF00001"]["type"] == "Family"
        assert result["PF00001"]["details"] == "7 transmembrane receptor"

    def test_load_rvdb(self, tmp_path):
        data = [{"RVDB_TargetID": "RVDB001", "RVDB_Description": "RNA polymerase"}]
        json_file = tmp_path / "rvdb.json"
        json_file.write_text(json.dumps(data))
        result = HmmMetadataLoader.load(str(json_file), "RVDB")
        assert "RVDB001" in result

    def test_load_vfam(self, tmp_path):
        data = [{"Vfam_TargetID": "Vfam001", "Vfam_Description": "Coat protein"}]
        json_file = tmp_path / "vfam.json"
        json_file.write_text(json.dumps(data))
        result = HmmMetadataLoader.load(str(json_file), "Vfam")
        assert "Vfam001" in result

    def test_unknown_db_returns_empty(self, tmp_path):
        json_file = tmp_path / "x.json"
        json_file.write_text("[]")
        result = HmmMetadataLoader.load(str(json_file), "Unknown")
        assert result == {}

    def test_missing_file_returns_empty(self):
        result = HmmMetadataLoader.load("/nonexistent/path.json", "EggNOG")
        assert result == {}

    def test_items_without_id_key_are_skipped(self, tmp_path):
        data = [
            {"EggNOG_Description": "no id here"},
            {"EggNOG_TargetID": "COG003", "EggNOG_Description": "has id"},
        ]
        json_file = tmp_path / "eggnog.json"
        json_file.write_text(json.dumps(data))
        result = HmmMetadataLoader.load(str(json_file), "EggNOG")
        assert len(result) == 1
        assert "COG003" in result


# ---------------------------------------------------------------------------
# HmmSequencePreparer
# ---------------------------------------------------------------------------

class TestHmmSequencePreparer:
    def test_empty_list_returns_none(self):
        block, orf_map = HmmSequencePreparer.prepare([])
        assert block is None
        assert orf_map == {}

    def test_nuc_seq_without_orfs_returns_none(self):
        block, orf_map = HmmSequencePreparer.prepare([make_nuc_seq()])
        assert block is None
        assert orf_map == {}

    def test_orf_with_empty_aa_sequence_is_skipped(self):
        seq = make_nuc_seq()
        seq.orfs.append(make_orf("orf1", aa_sequence=""))
        block, orf_map = HmmSequencePreparer.prepare([seq])
        assert block is None
        assert orf_map == {}

    def test_valid_orf_produces_seq_block(self):
        seq = make_nuc_seq()
        orf = make_orf("orf1", aa_sequence="MKVLS")
        seq.orfs.append(orf)
        block, orf_map = HmmSequencePreparer.prepare([seq])
        assert block is not None
        assert "orf1" in orf_map
        assert orf_map["orf1"] is orf

    def test_multiple_orfs_all_indexed(self):
        seq = make_nuc_seq()
        for i in range(3):
            seq.orfs.append(make_orf(f"orf{i}", aa_sequence="MKVLS"))
        block, orf_map = HmmSequencePreparer.prepare([seq])
        assert len(orf_map) == 3
        assert block is not None

    def test_mixed_orfs_only_valid_indexed(self):
        seq = make_nuc_seq()
        seq.orfs.append(make_orf("good", aa_sequence="MKVLS"))
        seq.orfs.append(make_orf("empty", aa_sequence=""))
        block, orf_map = HmmSequencePreparer.prepare([seq])
        assert "good" in orf_map
        assert "empty" not in orf_map


# ---------------------------------------------------------------------------
# HmmResultAttacher
# ---------------------------------------------------------------------------

class TestHmmResultAttacher:
    def _hit(self, query="HMM001", target="orf1", score=100.0,
             evalue=1e-10, env_from=1, env_to=50):
        return (query, target, score, evalue, env_from, env_to)

    def test_attach_adds_domain(self):
        orf = make_orf("orf1")
        metadata = {"HMM001": {"description": "Viral RdRp", "type": "Family", "details": "RNA dep"}}
        count = HmmResultAttacher.attach([self._hit()], {"orf1": orf}, metadata, "RVDB")
        assert count == 1
        assert len(orf.domains) == 1
        d = orf.domains[0]
        assert d.database == "RVDB"
        assert d.target == "HMM001"
        assert d.score == 100.0
        assert d.description == "Viral RdRp"
        assert d.type == "Family"
        assert d.details == "RNA dep"

    def test_domain_length_computed_correctly(self):
        orf = make_orf("orf1")
        HmmResultAttacher.attach([self._hit(env_from=10, env_to=60)], {"orf1": orf}, {}, "RVDB")
        assert orf.domains[0].length == 50

    def test_unknown_target_skipped(self):
        count = HmmResultAttacher.attach([self._hit(target="ghost")], {}, {}, "RVDB")
        assert count == 0

    def test_missing_metadata_yields_empty_strings(self):
        orf = make_orf("orf1")
        HmmResultAttacher.attach([self._hit()], {"orf1": orf}, {}, "RVDB")
        d = orf.domains[0]
        assert d.description == ""
        assert d.type == ""
        assert d.details == ""

    def test_multiple_hits_attached(self):
        orf = make_orf("orf1")
        hits = [self._hit(query="HMM001"), self._hit(query="HMM002")]
        count = HmmResultAttacher.attach(hits, {"orf1": orf}, {}, "RVDB")
        assert count == 2
        assert len(orf.domains) == 2


# ---------------------------------------------------------------------------
# HmmViralFlagSetter
# ---------------------------------------------------------------------------

class TestHmmViralFlagSetter:
    def test_rvdb_domain_flags_sequence(self):
        seq = make_nuc_seq()
        orf = make_orf("orf1")
        orf.domains.append(make_domain("RVDB"))
        seq.orfs.append(orf)
        count = HmmViralFlagSetter.flag([seq])
        assert seq.is_viral is True
        assert count == 1

    def test_vfam_domain_flags_sequence(self):
        seq = make_nuc_seq()
        orf = make_orf("orf1")
        orf.domains.append(make_domain("Vfam"))
        seq.orfs.append(orf)
        count = HmmViralFlagSetter.flag([seq])
        assert seq.is_viral is True
        assert count == 1

    def test_eggnon_domain_flags_sequence(self):
        seq = make_nuc_seq()
        orf = make_orf("orf1")
        orf.domains.append(make_domain("EggNOG"))
        seq.orfs.append(orf)
        count = HmmViralFlagSetter.flag([seq])
        assert seq.is_viral is True
        assert count == 1

    def test_pfam_domain_does_not_flag(self):
        seq = make_nuc_seq()
        orf = make_orf("orf1")
        orf.domains.append(make_domain("Pfam"))
        seq.orfs.append(orf)
        count = HmmViralFlagSetter.flag([seq])
        assert seq.is_viral is False
        assert count == 0

    def test_already_viral_not_counted(self):
        seq = make_nuc_seq()
        seq.is_viral = True
        orf = make_orf("orf1")
        orf.domains.append(make_domain("RVDB"))
        seq.orfs.append(orf)
        count = HmmViralFlagSetter.flag([seq])
        assert count == 0

    def test_no_domains_not_flagged(self):
        seq = make_nuc_seq()
        seq.orfs.append(make_orf("orf1"))
        count = HmmViralFlagSetter.flag([seq])
        assert seq.is_viral is False
        assert count == 0

    def test_empty_list(self):
        assert HmmViralFlagSetter.flag([]) == 0

    def test_multiple_sequences_partial_flag(self):
        seq1 = make_nuc_seq("seq1")
        orf1 = make_orf("orf1")
        orf1.domains.append(make_domain("RVDB"))
        seq1.orfs.append(orf1)

        seq2 = make_nuc_seq("seq2")
        seq2.orfs.append(make_orf("orf2"))

        count = HmmViralFlagSetter.flag([seq1, seq2])
        assert count == 1
        assert seq1.is_viral is True
        assert seq2.is_viral is False


# ---------------------------------------------------------------------------
# HmmSearcher
# ---------------------------------------------------------------------------

class TestHmmSearcher:
    def test_default_attributes(self):
        s = HmmSearcher()
        assert s.cpus == 0
        assert s.score_threshold == 50.0

    def test_custom_attributes(self):
        s = HmmSearcher(cpus=4, score_threshold=75.0)
        assert s.cpus == 4
        assert s.score_threshold == 75.0

    def _make_top_hits(self, query_name: bytes, hits):
        top_hits = MagicMock()
        top_hits.query.name = query_name
        top_hits.__iter__ = MagicMock(return_value=iter(hits))
        return top_hits

    def _make_hit(self, name: bytes, included: bool, domains):
        hit = MagicMock()
        hit.included = included
        hit.name = name
        hit.domains = domains
        return hit

    def _make_domain(self, score: float, env_from=1, env_to=50, i_evalue=1e-10):
        d = MagicMock()
        d.score = score
        d.i_evalue = i_evalue
        d.env_from = env_from
        d.env_to = env_to
        return d

    @patch("viralquest.hmm.pyhmmer.plan7.HMMFile")
    @patch("viralquest.hmm.pyhmmer.hmmsearch")
    def test_search_returns_hits_above_threshold(self, mock_hmmsearch, mock_hmmfile):
        domain = self._make_domain(score=80.0)
        hit = self._make_hit(b"orf1", included=True, domains=[domain])
        top_hits = self._make_top_hits(b"HMM001", [hit])
        mock_hmmsearch.return_value = [top_hits]

        results = HmmSearcher(score_threshold=50.0).search(MagicMock(), "/fake.hmm")
        assert results == [("HMM001", "orf1", 80.0, 1e-10, 1, 50)]

    @patch("viralquest.hmm.pyhmmer.plan7.HMMFile")
    @patch("viralquest.hmm.pyhmmer.hmmsearch")
    def test_search_filters_domains_below_threshold(self, mock_hmmsearch, mock_hmmfile):
        domain = self._make_domain(score=20.0)
        hit = self._make_hit(b"orf1", included=True, domains=[domain])
        top_hits = self._make_top_hits(b"HMM001", [hit])
        mock_hmmsearch.return_value = [top_hits]

        results = HmmSearcher(score_threshold=50.0).search(MagicMock(), "/fake.hmm")
        assert results == []

    @patch("viralquest.hmm.pyhmmer.plan7.HMMFile")
    @patch("viralquest.hmm.pyhmmer.hmmsearch")
    def test_search_skips_excluded_hits(self, mock_hmmsearch, mock_hmmfile):
        hit = self._make_hit(b"orf1", included=False, domains=[self._make_domain(80.0)])
        top_hits = self._make_top_hits(b"HMM001", [hit])
        mock_hmmsearch.return_value = [top_hits]

        results = HmmSearcher().search(MagicMock(), "/fake.hmm")
        assert results == []

    @patch("viralquest.hmm.pyhmmer.plan7.HMMFile")
    @patch("viralquest.hmm.pyhmmer.hmmsearch")
    def test_search_decodes_str_names(self, mock_hmmsearch, mock_hmmfile):
        domain = self._make_domain(score=80.0)
        hit = self._make_hit(b"orf1", included=True, domains=[domain])
        top_hits = self._make_top_hits(b"HMM001", [hit])
        mock_hmmsearch.return_value = [top_hits]

        results = HmmSearcher(score_threshold=50.0).search(MagicMock(), "/fake.hmm")
        query_name, target_name = results[0][0], results[0][1]
        assert isinstance(query_name, str)
        assert isinstance(target_name, str)

    def test_search_parallel_aggregates_results(self):
        searcher = HmmSearcher()
        fake_hit = ("HMM001", "orf1", 80.0, 1e-10, 1, 50)
        with patch.object(searcher, "search", return_value=[fake_hit]):
            results = searcher.search_parallel(MagicMock(), ["/a.hmm", "/b.hmm"])
        assert results["/a.hmm"] == [fake_hit]
        assert results["/b.hmm"] == [fake_hit]

    def test_search_parallel_handles_exception(self):
        searcher = HmmSearcher()
        with patch.object(searcher, "search", side_effect=RuntimeError("fail")):
            results = searcher.search_parallel(MagicMock(), ["/bad.hmm"])
        assert results["/bad.hmm"] == []
