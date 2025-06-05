import tarfile
import zipfile
import pytest
from nplinker.metabolomics.gnps import GNPSDownloader
from nplinker.metabolomics.gnps import GNPSFormat


def test_invalid_gnps_version(tmpdir):
    with pytest.raises(ValueError, match="Invalid GNPS version '3'"):
        GNPSDownloader("0ad6535e34d449788f297e712f43068a", tmpdir, "3")


def test_unknown_workflow(tmpdir):
    with pytest.raises(ValueError, match="Unknown workflow type for GNPS task .*"):
        GNPSDownloader("0ad6535e34d449788f297e712f43068a", tmpdir)


@pytest.mark.parametrize(
    "gnps_version, task_id, filename",
    [
        [
            "1",
            "92036537c21b44c29e509291e53f6382",
            GNPSFormat.FBMN.value + "-92036537c21b44c29e509291e53f6382.zip",
        ],
        [
            "1",
            "c22f44b14a3d450eb836d607cb9521bb",
            GNPSFormat.SNETS.value + "-c22f44b14a3d450eb836d607cb9521bb.zip",
        ],
        [
            "1",
            "189e8bf16af145758b0a900f1c44ff4a",
            GNPSFormat.SNETSV2.value + "-189e8bf16af145758b0a900f1c44ff4a.zip",
        ],
        [
            "2",
            "206a7b40b7ed41c1ae6b4fbd2def3636",
            "206a7b40b7ed41c1ae6b4fbd2def3636.tar",
        ],
        ["2", "2014f321d72542afb5216c932e0d5079", "2014f321d72542afb5216c932e0d5079.tar"],
    ],
)
def test_get_download_file(gnps_version, task_id, filename, tmpdir, gnps_website_is_down):
    if gnps_website_is_down:
        pytest.skip("GNPS website is down: https://gnps.ucsd.edu")
    downloader = GNPSDownloader(task_id, tmpdir, gnps_version)
    assert downloader.get_download_file() == tmpdir / filename


@pytest.mark.parametrize(
    "task_id",
    [
        "92036537c21b44c29e509291e53f6382",
        "c22f44b14a3d450eb836d607cb9521bb",
        "189e8bf16af145758b0a900f1c44ff4a",
    ],
)
def test_get_task_id(task_id, tmpdir, gnps_website_is_down):
    if gnps_website_is_down:
        pytest.skip("GNPS website is down: https://gnps.ucsd.edu")
    downloader = GNPSDownloader(task_id, tmpdir)
    assert downloader.get_task_id() == task_id


@pytest.mark.parametrize(
    "gnps_version, task_id, url",
    [
        [
            "1",
            "92036537c21b44c29e509291e53f6382",
            GNPSDownloader.GNPS_DATA_DOWNLOAD_URL_FBMN.format("92036537c21b44c29e509291e53f6382"),
        ],
        [
            "1",
            "c22f44b14a3d450eb836d607cb9521bb",
            GNPSDownloader.GNPS_DATA_DOWNLOAD_URL.format("c22f44b14a3d450eb836d607cb9521bb"),
        ],
        [
            "1",
            "189e8bf16af145758b0a900f1c44ff4a",
            GNPSDownloader.GNPS_DATA_DOWNLOAD_URL.format("189e8bf16af145758b0a900f1c44ff4a"),
        ],
        [
            "2",
            "206a7b40b7ed41c1ae6b4fbd2def3636",
            GNPSDownloader.GNPS2_DATA_DOWNLOAD_URL.format("206a7b40b7ed41c1ae6b4fbd2def3636"),
        ],
        [
            "2",
            "2014f321d72542afb5216c932e0d5079",
            GNPSDownloader.GNPS2_DATA_DOWNLOAD_URL.format("2014f321d72542afb5216c932e0d5079"),
        ],
    ],
)
def test_get_url(gnps_version, task_id, url, tmpdir, gnps_website_is_down):
    if gnps_website_is_down:
        pytest.skip("GNPS website is down: https://gnps.ucsd.edu")
    downloader = GNPSDownloader(task_id, tmpdir, gnps_version)
    assert downloader.get_url() == url


@pytest.mark.parametrize(
    "task_id, workflow",
    [
        ["92036537c21b44c29e509291e53f6382", GNPSFormat.FBMN],
        ["c22f44b14a3d450eb836d607cb9521bb", GNPSFormat.SNETS],
        ["189e8bf16af145758b0a900f1c44ff4a", GNPSFormat.SNETSV2],
    ],
)
def test_download_gnps1(task_id, workflow, tmpdir, gnps_zip_files, gnps_website_is_down):
    if gnps_website_is_down:
        pytest.skip("GNPS website is down: https://gnps.ucsd.edu")
    downloader = GNPSDownloader(task_id, tmpdir, gnps_version="1")
    downloader.download()
    actual = zipfile.ZipFile(downloader.get_download_file())
    actual_names = actual.namelist()
    expected = zipfile.ZipFile(gnps_zip_files[workflow])
    expected_names = expected.namelist()
    assert actual_names == expected_names


@pytest.mark.parametrize(
    "task_id, workflow",
    [
        ["2014f321d72542afb5216c932e0d5079", GNPSFormat.GNPS2FBMN],
        ["206a7b40b7ed41c1ae6b4fbd2def3636", GNPSFormat.GNPS2CN],
    ],
)
def test_download_gnps2(task_id, workflow, tmpdir, gnps2_tar_files, gnps2_website_is_down):
    if gnps2_website_is_down:
        pytest.skip("GNPS2 website is down: https://gnps2.org")
    downloader = GNPSDownloader(task_id, tmpdir, gnps_version="2")
    downloader.download()
    actual = tarfile.open(downloader.get_download_file())
    actual_names = actual.getnames()
    expected = tarfile.open(gnps2_tar_files[workflow])
    expected_names = expected.getnames()
    assert actual_names == expected_names
