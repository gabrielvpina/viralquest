import pytest
import subprocess
from viralquest.parser import FastaParser, Cap3Runner

def test_parser_read_fasta(tmp_path):
    
    mock_file = tmp_path / "test.fasta"
    mock_file.write_text(">seq1\nATGCGT\n>seq2\nCGTANC")

    parser = FastaParser(str(mock_file))
    parser.read_input_file()

    # assert checks
    assert len(parser.sequences) == 2
    assert parser.sequences[0].id == "seq1"
    assert parser.sequences[0].length == 6
    assert parser.sequences[1].n_count == 1


def test_parser_read_protein_fasta(tmp_path):

    mock_file_aa = tmp_path / "test_aa.fasta"
    mock_file_aa.write_text(">seq1\nMKWVT")

    parser = FastaParser(str(mock_file_aa))

    with pytest.raises(ValueError, match="invalid"):
        parser.read_input_file()



def test_cap3_runner_rejects_invalid_fasta(tmp_path):

    mock_file = tmp_path / "invalid_test.fasta"
    mock_file.write_text(">seq1\nATG1234@GTC")

    parser = FastaParser(str(mock_file))

    with pytest.raises(ValueError, match="invalid to the analysis"):
        parser.read_input_file()



def test_cap3_runner_reads_valid_fasta_successfully(tmp_path, monkeypatch):
    mock_file = tmp_path / "valid_test.fasta"
    mock_out = tmp_path / "results_cap3"
    mock_file.write_text(">seq1\nATGCGTACGTAC\n>seq2\nCGTACGTACGTA")

    runner = Cap3Runner(file_path=mock_file, outdir=mock_out)

    def mock_run(*args, **kwargs):

        log_file = runner.working_fasta.with_suffix(".cap3.log")
        log_file.write_text("CAP3 mock output")
        
        base = runner.working_fasta.name
        
        (runner.working_fasta.parent / f"{base}.cap.contigs").write_text(">link1\nAMOSTRA")
        (runner.working_fasta.parent / f"{base}.cap.singlets").write_text(">link2\nAMOSTRA")
        (runner.working_fasta.parent / f"{base}.cap.info").write_text("info")
        (runner.working_fasta.parent / f"{base}.cap.ace").write_text("ace")

        class MockResult:
            stdout = "CAP3 finished successfully"
        return MockResult()

    import subprocess
    monkeypatch.setattr(subprocess, "run", mock_run)

    resultado = runner.cap3_runner() 
    
    # Asserts
    assert runner.working_fasta.exists() is True 
    assert runner.working_fasta.parent == mock_out
    assert resultado.is_successful is True
