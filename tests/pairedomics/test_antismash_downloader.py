from pathlib import Path
from nplinker.pairedomics.podp_antismash_downloader import (
    GenomeStatus,
    podp_download_and_extract_antismash_data,
    _get_genome_status_log)
from nplinker.utils import list_files

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

def test_default(tmp_path):
    download_root = tmp_path / "download"
    extract_root = tmp_path / "extracted"
    genome_records = [
        {
            "genome_ID": {
                "genome_type": "genome",
                "JGI_Genome_ID": "2515154188",
                "RefSeq_accession": "GCF_000514875.1"
                },
            "genome_label": "Salinispora arenicola CNX508"},
        {
            "genome_ID": {
                "genome_type": "genome",
                "JGI_Genome_ID": "2515154177",
                "RefSeq_accession": "GCF_000514515.1"
                },
            "genome_label": "Salinispora pacifica CNT029"}
    ]
    genome_status_file = Path(download_root, "genome_status.csv")

    podp_download_and_extract_antismash_data(genome_records, download_root, extract_root)

    archive1 = download_root / "GCF_000514875.1.zip"
    extracted_folder1 = extract_root / "antismash" / "GCF_000514875.1"
    extracted_files1 = list_files(extracted_folder1, keep_parent=False)
    archive2 = download_root / "GCF_000514515.1.zip"
    extracted_folder2 = extract_root / "antismash" / "GCF_000514515.1"
    extracted_files2 = list_files(extracted_folder2, keep_parent=False)
    genome_status = _get_genome_status_log(genome_status_file)

    assert (archive1.exists() and archive2.exists())
    assert (archive1.is_file() and archive2.is_file())
    assert (extracted_folder1.exists() and extracted_folder2.exists())
    assert all(Path(extracted_folder1, extracted_file).is_file() for extracted_file in extracted_files1)
    assert all(Path(extracted_folder2, extracted_file).is_file() for extracted_file in extracted_files2)
    assert genome_status_file.is_file()
    assert len(genome_status) == 2
