import zipfile
import pytest
from nplinker.metabolomics.gnps import GNPSDownloader
from nplinker.metabolomics.gnps import GNPSFormat


@pytest.fixture(scope="module", autouse=True)
def setup_with_fixture(gnps_website_is_down):
    if gnps_website_is_down:
        pytest.skip(
            "GNPS website is down, skipping all tests in this module!", allow_module_level=True
        )


def test_unknown_workflow(tmpdir):
    with pytest.raises(ValueError, match="Unknown workflow type for GNPS task .*"):
        GNPSDownloader("0ad6535e34d449788f297e712f43068a", tmpdir)


@pytest.mark.parametrize(
    "task_id, expected",
    [
        ["92036537c21b44c29e509291e53f6382", GNPSFormat.FBMN],
        ["c22f44b14a3d450eb836d607cb9521bb", GNPSFormat.SNETS],
        ["189e8bf16af145758b0a900f1c44ff4a", GNPSFormat.SNETSV2],
    ],
)
def test_supported_workflows(task_id, expected, tmpdir):
    downloader = GNPSDownloader(task_id, tmpdir)
    assert downloader.gnps_format == expected


@pytest.mark.parametrize(
    "task_id, filename",
    [
        [
            "92036537c21b44c29e509291e53f6382",
            GNPSFormat.FBMN.value + "-92036537c21b44c29e509291e53f6382.zip",
        ],
        [
            "c22f44b14a3d450eb836d607cb9521bb",
            GNPSFormat.SNETS.value + "-c22f44b14a3d450eb836d607cb9521bb.zip",
        ],
        [
            "189e8bf16af145758b0a900f1c44ff4a",
            GNPSFormat.SNETSV2.value + "-189e8bf16af145758b0a900f1c44ff4a.zip",
        ],
    ],
)
def test_get_download_file(task_id, filename, tmpdir):
    downloader = GNPSDownloader(task_id, tmpdir)
    assert downloader.get_download_file() == tmpdir / filename


@pytest.mark.parametrize(
    "task_id",
    [
        "92036537c21b44c29e509291e53f6382",
        "c22f44b14a3d450eb836d607cb9521bb",
        "189e8bf16af145758b0a900f1c44ff4a",
    ],
)
def test_get_task_id(task_id, tmpdir):
    downloader = GNPSDownloader(task_id, tmpdir)
    assert downloader.get_task_id() == task_id


@pytest.mark.parametrize(
    "task_id, url",
    [
        [
            "92036537c21b44c29e509291e53f6382",
            GNPSDownloader.GNPS_DATA_DOWNLOAD_URL_FBMN.format("92036537c21b44c29e509291e53f6382"),
        ],
        [
            "c22f44b14a3d450eb836d607cb9521bb",
            GNPSDownloader.GNPS_DATA_DOWNLOAD_URL.format("c22f44b14a3d450eb836d607cb9521bb"),
        ],
        [
            "189e8bf16af145758b0a900f1c44ff4a",
            GNPSDownloader.GNPS_DATA_DOWNLOAD_URL.format("189e8bf16af145758b0a900f1c44ff4a"),
        ],
    ],
)
def test_get_url(task_id, url, tmpdir):
    downloader = GNPSDownloader(task_id, tmpdir)
    assert downloader.get_url() == url


@pytest.mark.parametrize(
    "task_id, workflow",
    [
        ["92036537c21b44c29e509291e53f6382", GNPSFormat.FBMN],
        ["c22f44b14a3d450eb836d607cb9521bb", GNPSFormat.SNETS],
        ["189e8bf16af145758b0a900f1c44ff4a", GNPSFormat.SNETSV2],
    ],
)
def test_downloads_file(task_id, workflow, tmpdir, gnps_zip_files):
    downloader = GNPSDownloader(task_id, tmpdir)
    downloader.download()
    actual = zipfile.ZipFile(downloader.get_download_file())
    actual_names = actual.namelist()
    expected = zipfile.ZipFile(gnps_zip_files[workflow])
    expected_names = [x.filename for x in expected.filelist if x.compress_size > 0]
    assert all(item in actual_names for item in expected_names)
