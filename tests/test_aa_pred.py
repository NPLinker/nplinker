import pytest
from Bio import SeqIO
from nplinker.genomics.aa_pred import AntiSmash5Record, predict_aa
from . import DATA_DIR


ANTISMASH_FILE = DATA_DIR / "antismash_v5_GCF_000016425.1_NC_009380.1.region017.gbk"

def test_predict_aa():
    pred = list(predict_aa(ANTISMASH_FILE))
    print(pred)
    assert len(pred) == 22
    assert ("ala", 0.0) in pred
    assert ("gly", 1.0) in pred
    assert ("val", 1.0) in pred


# Test class AntiSmash5Record
@pytest.fixture()
def antismash5_record():
        record = AntiSmash5Record(SeqIO.read(ANTISMASH_FILE, "genbank"))
        yield record

def test_get_prob(antismash5_record):
    assert antismash5_record.get_prob("ala") == 0.0
    assert antismash5_record.get_prob("gly") == 1.0
    assert antismash5_record.get_prob("val") == 1.0

def test_get_spec(antismash5_record):
    aa = list(antismash5_record.get_spec())
    assert len(aa) == 2
    assert 'gly' in aa
    assert 'val' in aa
