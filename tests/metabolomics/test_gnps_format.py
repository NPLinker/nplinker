import pytest
from nplinker.metabolomics.gnps import gnps_format_from_archive
from nplinker.metabolomics.gnps import gnps_format_from_file_mapping
from nplinker.metabolomics.gnps import gnps_format_from_task_id
from nplinker.metabolomics.gnps import GNPSFormat


@pytest.mark.parametrize(
    "task_id, expected",
    [
        ["92036537c21b44c29e509291e53f6382", GNPSFormat.FBMN],
        ["c22f44b14a3d450eb836d607cb9521bb", GNPSFormat.SNETS],
        ["189e8bf16af145758b0a900f1c44ff4a", GNPSFormat.SNETSV2],
        ["0ad6535e34d449788f297e712f43068a", GNPSFormat.Unknown],
    ],
)
def test_gnps_format_from_task_id(task_id: str, expected: GNPSFormat, gnps_website_is_down):
    if gnps_website_is_down:
        pytest.skip("GNPS website is down")
    actual = gnps_format_from_task_id(task_id)
    assert actual is expected


@pytest.mark.parametrize(
    "workflow", [GNPSFormat.FBMN, GNPSFormat.SNETS, GNPSFormat.SNETSV2, GNPSFormat.Unknown]
)
def test_gnps_format_from_archive(workflow: str, gnps_zip_files):
    actual = gnps_format_from_archive(gnps_zip_files[workflow])
    assert actual is workflow


@pytest.mark.parametrize("workflow", [GNPSFormat.FBMN, GNPSFormat.SNETS, GNPSFormat.SNETSV2])
def test_gnps_format_from_file_mapping(workflow: str, gnps_file_mappings_files):
    actual = gnps_format_from_file_mapping(gnps_file_mappings_files[workflow])
    assert actual is workflow
