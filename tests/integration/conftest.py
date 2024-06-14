from __future__ import annotations
import os
import shutil
import zipfile
from os import PathLike
from pathlib import Path
import httpx
import pytest
from rich.progress import Progress
from . import DATA_DIR


# The DOI of the test dataset to download, see https://zenodo.org/records/10822604
# NOTE: If you update the dataset doi, make sure remove the cached dataset "nplinker_local_mode_example.zip"
# from `DATA_DIR`. The code won't verify if the cached dataset is the same as the new dataset doi.
dataset_doi = "10.5281/zenodo.10822604"
dataset_url = (
    f"https://zenodo.org/records/{dataset_doi.split('.')[-1]}/files/nplinker_local_mode_example.zip"
)


@pytest.fixture(scope="module")
def root_dir(tmp_path_factory):
    """Set up the NPLinker root directory for the local mode example dataset."""
    temp_dir = tmp_path_factory.mktemp("nplinker_integration_test")
    nplinker_root_dir = temp_dir / "nplinker_local_mode_example"

    # Download the dataset and extract it
    if os.path.exists(nplinker_root_dir):
        shutil.rmtree(nplinker_root_dir)
    dataset = DATA_DIR / "nplinker_local_mode_example.zip"
    if not dataset.exists():
        download_archive(dataset_url, DATA_DIR)
    # the extracted directory is named "nplinker_local_mode_example"
    with zipfile.ZipFile(dataset, "r") as zip_ref:
        zip_ref.extractall(temp_dir)

    # Return the root directory
    yield str(nplinker_root_dir)

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

            with Progress() as progress:
                task = progress.add_task(f"[hot_pink]Downloading {fpath.name}", total=total)
                for chunk in response.iter_bytes():
                    fh.write(chunk)
                    progress.update(task, advance=len(chunk))
