import json
from pathlib import Path
import tempfile
import pytest
from nplinker.pairedomics import GENOME_STATUS_FILENAME
from nplinker.pairedomics import GenomeStatus
from nplinker.pairedomics import podp_download_and_extract_antismash_data
from nplinker.utils import list_files


@pytest.fixture
def download_root(tmp_path):
    return tmp_path / "download"


@pytest.fixture
def extract_root(tmp_path):
    return tmp_path / "extracted"


@pytest.fixture
def genome_status_file(download_root):
    return Path(download_root, GENOME_STATUS_FILENAME)


# Test `GenomeStatus` class
@pytest.mark.parametrize("params, expected", [
    (["genome1"], ["genome1", "", False, ""]),
    (["genome1", "refseq1", True, "/path/to/file"
      ], ["genome1", "refseq1", True, "/path/to/file"]),
])
def test_genome_status_init(params, expected):
    gs = GenomeStatus(*params)
    assert [
        gs.original_id, gs.resolved_refseq_id, gs.resolve_attempted,
        gs.bgc_path
    ] == expected


def test_genome_status_load_from_json(tmp_path):
    data = {
        "genome_status": [{
            "original_id": "genome1",
            "resolved_refseq_id": "refseq1",
            "resolve_attempted": True,
            "bgc_path": "/path/to/bgc1"
        }, {
            "original_id": "genome2",
            "resolved_refseq_id": "",
            "resolve_attempted": False,
            "bgc_path": ""
        }],
        "version":
        "1.0"
    }
    file_path = tmp_path / GENOME_STATUS_FILENAME
    with open(file_path, "w") as f:
        json.dump(data, f)
    genome_status_dict = GenomeStatus.load_from_json(file_path)

    assert len(genome_status_dict) == 2
    assert genome_status_dict["genome1"].original_id == "genome1"
    assert genome_status_dict["genome1"].resolved_refseq_id == "refseq1"
    assert genome_status_dict["genome1"].resolve_attempted is True
    assert genome_status_dict["genome1"].bgc_path == "/path/to/bgc1"
    assert genome_status_dict["genome2"].original_id == "genome2"
    assert genome_status_dict["genome2"].resolved_refseq_id == ""
    assert genome_status_dict["genome2"].resolve_attempted is False
    assert genome_status_dict["genome2"].bgc_path == ""


def test_genome_status_save_to_json(tmp_path):
    genome_status_dict = {
        "genome1": GenomeStatus("genome1", "refseq1", True, "/path/to/bgc1"),
        "genome2": GenomeStatus("genome2", "", False, "")
    }
    GenomeStatus.save_to_json(genome_status_dict, tmp_path)
    with open(tmp_path / GENOME_STATUS_FILENAME, "r") as f:
        loaded_data = json.load(f)

    assert loaded_data["version"] == "1.0"
    assert len(loaded_data["genome_status"]) == 2
    assert loaded_data["genome_status"][0]["original_id"] == "genome1"
    assert loaded_data["genome_status"][0]["resolved_refseq_id"] == "refseq1"
    assert loaded_data["genome_status"][0]["resolve_attempted"] is True
    assert loaded_data["genome_status"][0]["bgc_path"] == "/path/to/bgc1"
    assert loaded_data["genome_status"][1]["original_id"] == "genome2"
    assert loaded_data["genome_status"][1]["resolved_refseq_id"] == ""
    assert loaded_data["genome_status"][1]["resolve_attempted"] is False
    assert loaded_data["genome_status"][1]["bgc_path"] == ""


# Test `podp_download_and_extract_antismash_data` function
# with multiple records containing three types of genome IDs
def test_multiple_records(download_root, extract_root, genome_status_file):
    genome_records = [{
        "genome_ID": {
            "genome_type": "genome",
            "JGI_Genome_ID": "2515154188",
            "RefSeq_accession": "GCF_000514875.1",
            "GenBank_accession": "GCA_004799605.1"
        },
        "genome_label": "Salinispora arenicola CNX508"
    }, {
        "genome_ID": {
            "genome_type": "genome",
            "JGI_Genome_ID": "2515154177",
            "RefSeq_accession": "GCF_000514515.1",
            "GenBank_accession": "GCA_004799605.1"
        },
        "genome_label": "Salinispora pacifica CNT029"
    }]

    podp_download_and_extract_antismash_data(genome_records, download_root,
                                             extract_root)

    archive1 = download_root / "GCF_000514875.1.zip"
    extracted_folder1 = extract_root / "antismash" / "GCF_000514875.1"
    extracted_files1 = list_files(extracted_folder1, keep_parent=True)
    archive2 = download_root / "GCF_000514515.1.zip"
    extracted_folder2 = extract_root / "antismash" / "GCF_000514515.1"
    extracted_files2 = list_files(extracted_folder2, keep_parent=True)
    genome_status = GenomeStatus.load_from_json(genome_status_file)

    assert archive1.exists()
    assert archive2.exists()
    assert archive1.is_file()
    assert archive2.is_file()
    assert extracted_folder1.exists()
    assert extracted_folder2.exists()
    assert all(
        Path(extracted_file).is_file() for extracted_file in extracted_files1)
    assert all(
        Path(extracted_file).is_file() for extracted_file in extracted_files2)
    assert genome_status_file.is_file()
    assert len(genome_status) == 2


# Test `podp_download_and_extract_antismash_data` function
# when a genome record has no accession IDs
def test_missing_id(download_root, extract_root, genome_status_file):
    genome_records = [{
        "genome_ID": {
            "genome_type": "genome",
            "JGI_Genome_ID": "",
            "RefSeq_accession": ""
        },
        "genome_label": "Salinispora arenicola CNX508"
    }, {
        "genome_ID": {
            "genome_type": "genome",
            "JGI_Genome_ID": "2515154177",
            "RefSeq_accession": "GCF_000514515.1"
        },
        "genome_label": "Salinispora pacifica CNT029"
    }]

    podp_download_and_extract_antismash_data(genome_records, download_root,
                                             extract_root)

    archive1 = download_root / "GCF_000514875.1.zip"
    extracted_folder1 = extract_root / "antismash" / "GCF_000514875.1"
    archive2 = download_root / "GCF_000514515.1.zip"
    extracted_folder2 = extract_root / "antismash" / "GCF_000514515.1"
    genome_status = GenomeStatus.load_from_json(genome_status_file)

    assert (not archive1.exists() and archive2.exists())
    assert (not archive1.is_file() and archive2.is_file())
    assert (not extracted_folder1.exists() and extracted_folder2.exists())
    assert genome_status_file.is_file()
    assert len(genome_status) == 1


# Test `podp_download_and_extract_antismash_data` function
# when a genome record has already been downloaded and extracted
def test_caching(download_root, extract_root, genome_status_file, caplog):
    genome_records = [{
        "genome_ID": {
            "genome_type": "genome",
            "JGI_Genome_ID": "2515154188",
            "RefSeq_accession": "GCF_000514875.1"
        },
        "genome_label": "Salinispora arenicola CNX508"
    }]

    podp_download_and_extract_antismash_data(genome_records, download_root,
                                             extract_root)
    genome_status_old = GenomeStatus.load_from_json(genome_status_file)
    genome_obj = genome_status_old["GCF_000514875.1"]
    assert Path(genome_obj.bgc_path).exists()
    assert genome_obj.resolve_attempted
    podp_download_and_extract_antismash_data(genome_records, download_root,
                                             extract_root)
    assert f'Genome ID {genome_obj.original_id} already downloaded to {genome_obj.bgc_path}' in caplog.text
    assert f'Genome ID {genome_obj.original_id} skipped due to previous failure' not in caplog.text
    genome_status_new = GenomeStatus.load_from_json(genome_status_file)
    assert len(genome_status_old) == len(genome_status_new)


# Test `podp_download_and_extract_antismash_data` function
# when a genome record has a non existing accession ID
def test_failed_lookup(download_root, extract_root, genome_status_file):
    genome_records = [{
        "genome_ID": {
            "genome_type": "genome",
            "JGI_Genome_ID": "non_existing_ID"
        },
        "genome_label": "Salinispora arenicola CNX508"
    }]

    podp_download_and_extract_antismash_data(genome_records, download_root,
                                             extract_root)
    genome_status = GenomeStatus.load_from_json(genome_status_file)
    assert len(genome_status["non_existing_ID"].bgc_path) == 0
    assert genome_status["non_existing_ID"].resolve_attempted
    assert not (download_root / "non_existing_ID.zip").exists()
    assert not (extract_root / "antismash" / "non_existing_ID.zip").exists()


# Test `podp_download_and_extract_antismash_data` function
# when a genome record has only RefSeq accession ID
def test_refseq_id(download_root, extract_root, genome_status_file):
    genome_records = [{
        "genome_ID": {
            "genome_type": "genome",
            "RefSeq_accession": "GCF_000514875.1"
        },
        "genome_label": "Salinispora arenicola CNX508"
    }]

    podp_download_and_extract_antismash_data(genome_records, download_root,
                                             extract_root)

    genome_status = GenomeStatus.load_from_json(genome_status_file)
    genome_obj = genome_status["GCF_000514875.1"]
    archive = download_root / Path(str(genome_obj.resolved_refseq_id) + ".zip")
    extracted_folder = extract_root / "antismash" / genome_obj.resolved_refseq_id
    extracted_files = list_files(extracted_folder, keep_parent=False)

    assert archive.exists()
    assert archive.is_file()
    assert extracted_folder.exists()
    assert all(
        Path(extracted_folder, extracted_file).is_file()
        for extracted_file in extracted_files)
    assert genome_status_file.is_file()
    assert len(genome_status) == 1


# Test `podp_download_and_extract_antismash_data` function
# when a genome record has only GenBank accession ID
def test_genbank_id(download_root, extract_root, genome_status_file):
    genome_records = [{
        "genome_ID": {
            "genome_type": "genome",
            "GenBank_accession": "GCA_004799605.1"
        },
        "genome_label": "Halobacterium salinarum"
    }]

    podp_download_and_extract_antismash_data(genome_records, download_root,
                                             extract_root)

    genome_status = GenomeStatus.load_from_json(genome_status_file)
    genome_obj = genome_status["GCA_004799605.1"]
    archive = download_root / Path(str(genome_obj.resolved_refseq_id) + ".zip")
    extracted_folder = extract_root / "antismash" / genome_obj.resolved_refseq_id
    extracted_files = list_files(extracted_folder, keep_parent=False)

    assert archive.exists()
    assert archive.is_file()
    assert extracted_folder.exists()
    assert all(
        Path(extracted_folder, extracted_file).is_file()
        for extracted_file in extracted_files)
    assert genome_status_file.is_file()
    assert len(genome_status) == 1


# Test `podp_download_and_extract_antismash_data` function
# when a genome record has only JGI accession ID
def test_jgi_id(download_root, extract_root, genome_status_file):
    genome_records = [{
        "genome_ID": {
            "genome_type": "genome",
            "JGI_Genome_ID": "2506783052"
        },
        "genome_label": "Halophilic archaeon DL31"
    }]

    podp_download_and_extract_antismash_data(genome_records, download_root,
                                             extract_root)

    genome_status = GenomeStatus.load_from_json(genome_status_file)
    genome_obj = genome_status["2506783052"]
    archive = download_root / Path(str(genome_obj.resolved_refseq_id) + ".zip")
    extracted_folder = extract_root / "antismash" / genome_obj.resolved_refseq_id
    extracted_files = list_files(extracted_folder, keep_parent=False)

    assert archive.exists()
    assert archive.is_file()
    assert extracted_folder.exists()
    assert all(
        Path(extracted_folder, extracted_file).is_file()
        for extracted_file in extracted_files)
    assert genome_status_file.is_file()
    assert len(genome_status) == 1


# Test `podp_download_and_extract_antismash_data` function
# when a genome record has an existing RefSeq and JGI accession ID;
# verify that RefSeq is used
def test_refseq_jgi_id(download_root, extract_root, genome_status_file):
    genome_records = [{
        "genome_ID": {
            "genome_type": "genome",
            "RefSeq_accession": "GCF_000514875.1",
            "JGI_Genome_ID": "2506783052"
        }
    }]

    podp_download_and_extract_antismash_data(genome_records, download_root,
                                             extract_root)

    genome_status = GenomeStatus.load_from_json(genome_status_file)
    genome_obj = genome_status["GCF_000514875.1"]
    archive = download_root / Path(str(genome_obj.resolved_refseq_id) + ".zip")
    extracted_folder = extract_root / "antismash" / genome_obj.resolved_refseq_id
    extracted_files = list_files(extracted_folder, keep_parent=False)

    assert archive.exists()
    assert archive.is_file()
    assert extracted_folder.exists()
    assert all(
        Path(extracted_folder, extracted_file).is_file()
        for extracted_file in extracted_files)
    assert genome_status_file.is_file()
    assert len(genome_status) == 1


# Test `podp_download_and_extract_antismash_data` function
# when a genome record has an existing RefSeq and GenBank accession ID;
# verify that RefSeq is used
def test_refseq_genbank_id(download_root, extract_root, genome_status_file):
    genome_records = [{
        "genome_ID": {
            "genome_type": "genome",
            "RefSeq_accession": "GCF_000514875.1",
            "GenBank_accession": "GCA_004799605.1"
        }
    }]

    podp_download_and_extract_antismash_data(genome_records, download_root,
                                             extract_root)

    genome_status = GenomeStatus.load_from_json(genome_status_file)
    genome_obj = genome_status["GCF_000514875.1"]
    archive = download_root / Path(str(genome_obj.resolved_refseq_id) + ".zip")
    extracted_folder = extract_root / "antismash" / genome_obj.resolved_refseq_id
    extracted_files = list_files(extracted_folder, keep_parent=False)

    assert archive.exists()
    assert archive.is_file()
    assert extracted_folder.exists()
    assert all(
        Path(extracted_folder, extracted_file).is_file()
        for extracted_file in extracted_files)
    assert genome_status_file.is_file()
    assert len(genome_status) == 1


# Test `podp_download_and_extract_antismash_data` function
# when a genome record has an existing GenBank and JGI accession ID;
# verify that GenBank is used
def test_genbank_jgi_id(download_root, extract_root, genome_status_file):
    genome_records = [{
        "genome_ID": {
            "genome_type": "genome",
            "GenBank_accession": "GCA_004799605.1",
            "JGI_Genome_ID": "2506783052"
        }
    }]

    podp_download_and_extract_antismash_data(genome_records, download_root,
                                             extract_root)

    genome_status = GenomeStatus.load_from_json(genome_status_file)
    genome_obj = genome_status["GCA_004799605.1"]
    archive = download_root / Path(str(genome_obj.resolved_refseq_id) + ".zip")
    extracted_folder = extract_root / "antismash" / genome_obj.resolved_refseq_id
    extracted_files = list_files(extracted_folder, keep_parent=False)

    assert archive.exists()
    assert archive.is_file()
    assert extracted_folder.exists()
    assert all(
        Path(extracted_folder, extracted_file).is_file()
        for extracted_file in extracted_files)
    assert genome_status_file.is_file()
    assert len(genome_status) == 1
