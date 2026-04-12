import pytest
from viralquest.parser import FastaParser

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
