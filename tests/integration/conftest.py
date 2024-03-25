from __future__ import annotations
import os
import shutil
import tempfile
import zipfile
from os import PathLike
from pathlib import Path
import httpx
from tqdm import tqdm
from . import DATA_DIR


# The DOI of the test dataset to download, see https://zenodo.org/records/10822604
# NOTE: If you update the dataset doi, make sure remove the cached dataset "nplinker_local_mode_example.zip"
# from `DATA_DIR`. The code won't verify if the cached dataset is the same as the new dataset doi.
dataset_doi = "10.5281/zenodo.10822604"
dataset_url = (
    f"https://zenodo.org/records/{dataset_doi.split('.')[-1]}/files/nplinker_local_mode_example.zip"
)

# The temporary directory for the test session
temp_dir = tempfile.gettempdir()
nplinker_root_dir = os.path.join(temp_dir, "nplinker_local_mode_example")


def pytest_sessionstart(session):
    """Pytest hook to run before the entire test session starts.

    This hook makes sure the temporary directory `nplinker_root_dir` is created before any test
    starts. When running tests in parallel, the creation operation is done by the master process,
    and worker processes are not allowed to do it.

    For more about this hook, see:
    1. https://docs.pytest.org/en/stable/reference.html#_pytest.hookspec.pytest_sessionstart
    2. https://github.com/pytest-dev/pytest-xdist/issues/271#issuecomment-826396320
    """
    workerinput = getattr(session.config, "workerinput", None)
    # It's master process or not running in parallell when `workerinput` is None.
    if workerinput is None:
        if os.path.exists(nplinker_root_dir):
            shutil.rmtree(nplinker_root_dir)
        dataset = DATA_DIR / "nplinker_local_mode_example.zip"
        if not dataset.exists():
            download_archive(dataset_url, DATA_DIR)
        with zipfile.ZipFile(dataset, "r") as zip_ref:
            zip_ref.extractall(temp_dir)
    # NPLinker setting `root_dir` must be a path that exists, so setting it to a temporary directory.
    os.environ["NPLINKER_ROOT_DIR"] = nplinker_root_dir
    # # Specify the config file via environment variable before importing nplinker in any test.
    os.environ["NPLINKER_CONFIG_FILE"] = str(DATA_DIR / "nplinker_local_mode.toml")


def pytest_sessionfinish(session):
    """Pytest hook to run after the entire test session finishes.

    This hook makes sure that temporary directory `nplinker_root_dir` is only removed after all
    tests finish. When running tests in parallel, the deletion operation is done by the master
    process, and worker processes are not allowed to do it.
    """
    workerinput = getattr(session.config, "workerinput", None)
    if workerinput is None:
        shutil.rmtree(nplinker_root_dir)


def download_archive(
    url: str,
    root: str | PathLike,
    http_method: str = "GET",
    allow_http_redirect: bool = True,
) -> None:
    """Download an archive file from the given URL.

    The output file name is determined by the URL and saved to the given root directory.

    Args:
        url: The URL of the file to download.
        root: The directory to save the file to.
        http_method: The HTTP method to use. Defaults to "GET".
        allow_http_redirect: Whether to allow HTTP redirects. Defaults to True.
    """
    fpath = Path(root) / Path(url).name

    # download the file
    with open(fpath, "wb") as fh:
        with httpx.stream(http_method, url, follow_redirects=allow_http_redirect) as response:
            response.raise_for_status()
            print(f"Downloading test dataset {url} to {root}")
            total = int(response.headers.get("Content-Length", 0))
            with tqdm(total=total, unit_scale=True, unit_divisor=1024, unit="B") as progress:
                num_bytes_downloaded = response.num_bytes_downloaded
                for chunk in response.iter_bytes():
                    fh.write(chunk)
                    progress.update(response.num_bytes_downloaded - num_bytes_downloaded)
                    num_bytes_downloaded = response.num_bytes_downloaded
