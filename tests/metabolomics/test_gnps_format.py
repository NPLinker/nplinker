import zipfile
import pytest

from nplinker.metabolomics.gnps.gnps_format import GNPSFormat
from nplinker.metabolomics.gnps.gnps_format import gnps_format_from_file_mapping
from nplinker.metabolomics.gnps.gnps_format import gnps_format_from_task_id
from nplinker.metabolomics.gnps.gnps_format import gnps_format_from_archive
from .. import DATA_DIR


@pytest.mark.parametrize("filename, expected", [
    [DATA_DIR / "nodes.tsv", GNPSFormat.AllFiles],
    [DATA_DIR / "nodes_mwe.csv", GNPSFormat.AllFiles],
    [DATA_DIR / "nodes_fbmn.csv", GNPSFormat.FBMN]
])
def test_identify_gnps_format(filename, expected):
    actual = gnps_format_from_file_mapping(filename, False)

    assert actual is expected


@pytest.mark.parametrize("task_id, expected", [
    ["92036537c21b44c29e509291e53f6382", GNPSFormat.FBMN],
    ["c22f44b14a3d450eb836d607cb9521bb", GNPSFormat.AllFiles]
])
def test_gnps_format_from_task_id(task_id: str, expected: GNPSFormat):
    actual = gnps_format_from_task_id(task_id)
    assert actual is expected

@pytest.mark.parametrize("archive_path, expected", [
    ["ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-92036537-download_cytoscape_data.zip", GNPSFormat.FBMN],
    ["ProteoSAFe-METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra.zip", GNPSFormat.AllFiles]
])
def test_gnps_format_from_archive(archive_path: str, expected: GNPSFormat):
    archive = zipfile.ZipFile(DATA_DIR / archive_path)

    actual = gnps_format_from_archive(archive)
    assert actual is expected
