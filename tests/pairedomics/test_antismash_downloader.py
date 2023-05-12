from pathlib import Path
from nplinker.pairedomics.podp_antismash_downloader import GenomeStatus, _get_genome_status_log


def test_GenomeStatus(tmp_path):
    genome_status_file = Path(tmp_path, "genome_status.csv")
    raw_genome_id1 = "GCF_000515175.1"
    raw_genome_id2 = "GCF_000514635.1"
    genome_obj1 = GenomeStatus(raw_genome_id1, "None")
    genome_obj2 = GenomeStatus(raw_genome_id2, "None")
    assert genome_obj1.resolved_refseq_id == ""
    assert genome_obj2.resolved_refseq_id == ""
    genome_status = _get_genome_status_log(genome_status_file)
    assert isinstance(genome_status, dict)
    assert len(genome_status) == 0
    genome_obj1.to_csv(genome_status_file)
    genome_obj2.to_csv(genome_status_file)
    genome_status = _get_genome_status_log(genome_status_file)
    assert isinstance(genome_status[raw_genome_id1], GenomeStatus)
    assert isinstance(genome_status[raw_genome_id2], GenomeStatus)
    assert genome_status[raw_genome_id1].original_id == raw_genome_id1
    assert genome_status[raw_genome_id2].original_id == raw_genome_id2