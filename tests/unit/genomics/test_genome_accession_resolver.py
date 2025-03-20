from unittest.mock import patch
import pytest
from nplinker.genomics.antismash.genome_accession_resolver import get_latest_assembly_accession
from nplinker.genomics.antismash.genome_accession_resolver import resolve_genome_accession
from nplinker.genomics.antismash.genome_accession_resolver import validate_assembly_acc


USER_AGENT = "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:86.0) Gecko/20100101 Firefox/86.0"


@pytest.fixture
def mock_resolvers():
    with (
        patch(
            "nplinker.genomics.antismash.genome_accession_resolver._resolve_refseq"
        ) as mock_refseq,
        patch(
            "nplinker.genomics.antismash.genome_accession_resolver._resolve_genbank"
        ) as mock_genbank,
        patch("nplinker.genomics.antismash.genome_accession_resolver._resolve_jgi") as mock_jgi,
    ):
        yield mock_refseq, mock_genbank, mock_jgi


# test get_latest_assembly_accession


def test_get_latest_assembly_accession_refseq_success():
    with patch(
        "nplinker.genomics.antismash.genome_accession_resolver._get_revision_history"
    ) as mock_get_revision_history:
        mock_get_revision_history.return_value = {
            "assembly_revisions": [
                {"refseq_accession": "GCF_000001.2", "release_date": "2023-01-01"},
                {"refseq_accession": "GCF_000001.1", "release_date": "2022-01-01"},
            ]
        }

        result = get_latest_assembly_accession("GCF_000001.1")

        assert result == "GCF_000001.2"
        mock_get_revision_history.assert_called_once_with("GCF_000001.1")


def test_get_latest_assembly_accession_genbank_success():
    with patch(
        "nplinker.genomics.antismash.genome_accession_resolver._get_revision_history"
    ) as mock_get_revision_history:
        mock_get_revision_history.return_value = {
            "assembly_revisions": [
                {"genbank_accession": "GCA_000002.2", "release_date": "2023-01-01"},
                {"genbank_accession": "GCA_000002.1", "release_date": "2022-01-01"},
            ]
        }

        result = get_latest_assembly_accession("GCA_000002.1")

        assert result == "GCA_000002.2"
        mock_get_revision_history.assert_called_once_with("GCA_000002.1")


def test_get_latest_assembly_accession_refseq_and_genbank():
    with patch(
        "nplinker.genomics.antismash.genome_accession_resolver._get_revision_history"
    ) as mock_get_revision_history:
        mock_get_revision_history.return_value = {
            "assembly_revisions": [
                {"refseq_accession": "GCF_000001.2", "release_date": "2023-01-01"},
                {"genbank_accession": "GCA_000002.2", "release_date": "2022-01-01"},
            ]
        }

        result = get_latest_assembly_accession("GCF_000001.1")

        assert result == "GCF_000001.2"
        mock_get_revision_history.assert_called_once_with("GCF_000001.1")


def test_get_latest_assembly_accession_no_valid_accession():
    with patch(
        "nplinker.genomics.antismash.genome_accession_resolver._get_revision_history"
    ) as mock_get_revision_history:
        mock_get_revision_history.return_value = {"assembly_revisions": []}

        with pytest.raises(ValueError, match="No valid genome accession found"):
            get_latest_assembly_accession("GCF_000001.1")

        mock_get_revision_history.assert_called_once_with("GCF_000001.1")


# test resolve_for_genome_accession


def test_resolve_genome_accession_refseq_success(mock_resolvers):
    mock_refseq, mock_genbank, mock_jgi = mock_resolvers
    mock_refseq.return_value = "GCF_000001.2"

    genome_id_data = {"RefSeq_accession": "GCF_000001.1"}
    result = resolve_genome_accession(genome_id_data)

    assert result == "GCF_000001.2"
    mock_refseq.assert_called_once_with("GCF_000001.1")
    mock_genbank.assert_not_called()
    mock_jgi.assert_not_called()


def test_resolve_genome_accession_genbank_success(mock_resolvers):
    mock_refseq, mock_genbank, mock_jgi = mock_resolvers
    mock_refseq.side_effect = Exception("RefSeq resolution failed")
    mock_genbank.return_value = "GCA_000002.2"

    genome_id_data = {
        "RefSeq_accession": "GCF_000001.1",
        "GenBank_accession": "GCA_000002.1",
    }
    result = resolve_genome_accession(genome_id_data)

    assert result == "GCA_000002.2"
    mock_refseq.assert_called_once_with("GCF_000001.1")
    mock_genbank.assert_called_once_with("GCA_000002.1")
    mock_jgi.assert_not_called()


def test_resolve_genome_accession_jgi_success(mock_resolvers):
    mock_refseq, mock_genbank, mock_jgi = mock_resolvers
    mock_refseq.side_effect = Exception("RefSeq resolution failed")
    mock_genbank.side_effect = Exception("GenBank resolution failed")
    mock_jgi.return_value = "GCF_000001.2"

    genome_id_data = {
        "RefSeq_accession": "GCF_000001.1",
        "GenBank_accession": "GCA_000002.1",
        "JGI_Genome_ID": "12345",
    }
    result = resolve_genome_accession(genome_id_data)

    assert result == "GCF_000001.2"
    mock_refseq.assert_called_once_with("GCF_000001.1")
    mock_genbank.assert_called_once_with("GCA_000002.1")
    mock_jgi.assert_called_once_with("12345")


def test_resolve_genome_accession_no_valid_accession(mock_resolvers):
    mock_refseq, mock_genbank, mock_jgi = mock_resolvers
    mock_refseq.side_effect = Exception("RefSeq resolution failed")
    mock_genbank.side_effect = Exception("GenBank resolution failed")
    mock_jgi.side_effect = Exception("JGI resolution failed")

    genome_id_data = {
        "RefSeq_accession": "GCF_000001.1",
        "GenBank_accession": "GCA_000002.1",
        "JGI_Genome_ID": "12345",
    }

    with pytest.raises(RuntimeError, match="No valid assembly accessions found"):
        resolve_genome_accession(genome_id_data)

    mock_refseq.assert_called_once_with("GCF_000001.1")
    mock_genbank.assert_called_once_with("GCA_000002.1")
    mock_jgi.assert_called_once_with("12345")


def test_resolve_genome_accession_missing_keys(mock_resolvers):
    mock_refseq, mock_genbank, mock_jgi = mock_resolvers

    genome_id_data = {}
    with pytest.raises(RuntimeError, match="No valid assembly accessions found"):
        resolve_genome_accession(genome_id_data)

    mock_refseq.assert_not_called()
    mock_genbank.assert_not_called()
    mock_jgi.assert_not_called()


def test_validate_assembly_acc_valid():
    # Test valid GenBank accession
    validate_assembly_acc("GCA_000002.1", "GenBank")
    # Test valid RefSeq accession
    validate_assembly_acc("GCF_000001.1", "RefSeq")


def test_validate_assembly_acc_invalid_type():
    # Test invalid accession type
    with pytest.raises(ValueError, match="Invalid genome assembly accession type 'InvalidType'"):
        validate_assembly_acc("GCF_000001.1", "InvalidType")


def test_validate_assembly_acc_invalid_prefix():
    # Test invalid prefix for GenBank accession
    with pytest.raises(ValueError, match="Invalid GenBank assembly accession"):
        validate_assembly_acc("GCF_000002.1", "GenBank")
    # Test invalid prefix for RefSeq accession
    with pytest.raises(ValueError, match="Invalid RefSeq assembly accession"):
        validate_assembly_acc("GCA_000001.1", "RefSeq")


def test_validate_assembly_acc_missing_version():
    # Test missing version number for GenBank accession
    with pytest.raises(ValueError, match="Invalid assembly accession \\(missing version number\\)"):
        validate_assembly_acc("GCA_000002", "GenBank")
    # Test missing version number for RefSeq accession
    with pytest.raises(ValueError, match="Invalid assembly accession \\(missing version number\\)"):
        validate_assembly_acc("GCF_000001", "RefSeq")
