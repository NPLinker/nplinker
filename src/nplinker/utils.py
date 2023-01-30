# Copyright 2021 The NPLinker Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import annotations
import csv
import os
from os import PathLike
import os.path
import bz2
import gzip
import hashlib
import lzma
import math
import sys
import tarfile
import urllib
import urllib.error
import urllib.request
import zipfile
from pathlib import Path
from typing import IO
from typing import Callable
from typing import Iterator
from tqdm import tqdm


# CG: it's only used by metabolomics, should move it there
# code to normalise peaks for spectral matching ("rosetta stone" stuff)
def sqrt_normalise(peaks):
    temp = []
    total = 0.0
    for mz, intensity in peaks:
        temp.append((mz, math.sqrt(intensity)))
        total += intensity
    norm_facc = math.sqrt(total)
    normalised_peaks = []
    for mz, intensity in temp:
        normalised_peaks.append((mz, intensity / norm_facc))
    return normalised_peaks


def find_delimiter(filename: str | PathLike) -> str:
    """Detect the delimiter for the given tabular file.

    Args:
        filename(str | PathLike): Path to tabular file.

    Returns:
        str: Detected delimiter character.

    Examples:
        >>> delim = find_delimiter("~/table.csv")
    """
    sniffer = csv.Sniffer()
    with open(filename, mode='rt', encoding='utf-8') as fp:
        delimiter = sniffer.sniff(fp.read(5000)).delimiter
    return delimiter
    
def get_headers(filename: str | PathLike) -> list[str]:
    """Function to read headers from tab or comma separated table.

    Args:
        filename(str): Path to the file to read the header from.

    Returns:
        list[str]: Columns names in header.
    """
    with open(os.fspath(filename)) as f:
        headers: str = f.readline().strip()
        dl: str = find_delimiter(filename)
        return headers.split(dl)

# Functions below are adapted from torchvision library,
# see: https://github.com/pytorch/vision/blob/main/torchvision/datasets/utils.py.
#
# BSD 3-Clause License
# Copyright (c) Soumith Chintala 2016
# You may obtain a copy of the License at
#    https://github.com/pytorch/vision/blob/main/LICENSE

USER_AGENT = "NPLinker"
# USER_AGENT = 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:86.0) Gecko/20100101 Firefox/86.0'


def _save_response_content(content: Iterator[bytes],
                           destination: str | Path,
                           length: int | None = None) -> None:
    with open(destination, "wb") as fh, tqdm(total=length) as pbar:
        for chunk in content:
            # filter out keep-alive new chunks
            if not chunk:
                continue

            fh.write(chunk)
            pbar.update(len(chunk))


def _urlretrieve(url: str,
                 filename: str | Path,
                 chunk_size: int = 1024 * 32) -> None:
    with urllib.request.urlopen(
            urllib.request.Request(url, headers={"User-Agent":
                                                 USER_AGENT})) as response:
        _save_response_content(iter(lambda: bytes(response.read(chunk_size)),
                                    b""),
                               filename,
                               length=response.length)


def calculate_md5(fpath: str | Path, chunk_size: int = 1024 * 1024) -> str:
    # Setting the `usedforsecurity` flag does not change anything about the functionality, but indicates that we are
    # not using the MD5 checksum for cryptography. This enables its usage in restricted environments like FIPS. Without
    # it torchvision.datasets is unusable in these environments since we perform a MD5 check everywhere.
    if sys.version_info >= (3, 9):
        md5 = hashlib.md5(usedforsecurity=False)
    else:
        md5 = hashlib.md5()
    with open(fpath, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            md5.update(chunk)
    return md5.hexdigest()


def check_md5(fpath: str | Path, md5: str) -> bool:
    return md5 == calculate_md5(fpath)


def check_integrity(fpath: str | Path, md5: str | None = None) -> bool:
    if not os.path.isfile(fpath):
        return False
    if md5 is None:
        return True
    return check_md5(fpath, md5)


def _get_redirect_url(url: str, max_hops: int = 3) -> str:
    initial_url = url
    headers = {"Method": "HEAD", "User-Agent": USER_AGENT}

    for _ in range(max_hops + 1):
        with urllib.request.urlopen(
                urllib.request.Request(url, headers=headers)) as response:
            if response.url == url or response.url is None:
                return url
            url = response.url

    raise RecursionError(
        f"Request to {initial_url} exceeded {max_hops} redirects. The last redirect points to {url}."
    )


def download_url(url: str,
                 root: str | Path,
                 filename: str | None = None,
                 md5: str | None = None,
                 max_redirect_hops: int = 3) -> None:
    """Download a file from a url and place it in root.

    Args:
        url (str): URL to download file from
        root (str): Directory to place downloaded file in
        filename (str, optional): Name to save the file under. If None, use the basename of the URL
        md5 (str, optional): MD5 checksum of the download. If None, do not check
        max_redirect_hops (int, optional): Maximum number of redirect hops allowed
    """
    root = Path(root).expanduser()
    if not filename:
        filename = Path(url).name
    fpath = root / filename

    Path.mkdir(root, exist_ok=True)

    # check if file is already present locally
    if check_integrity(fpath, md5):
        print("Using downloaded and verified file: " + str(fpath))
        return

    # expand redirect chain if needed
    url = _get_redirect_url(url, max_hops=max_redirect_hops)

    # download the file
    try:
        print("Downloading " + url + " to " + str(fpath))
        _urlretrieve(url, fpath)
    except (urllib.error.URLError, OSError) as e:  # type: ignore[attr-defined]
        if url[:5] == "https":
            url = url.replace("https:", "http:")
            print(
                "Failed download. Trying https -> http instead. Downloading " +
                url + " to " + str(fpath))
            _urlretrieve(url, fpath)
        else:
            raise e

    # check integrity of downloaded file
    if not check_integrity(fpath, md5):
        raise RuntimeError("File not found or corrupted, or md5 validation failed.")


def list_dirs(root: str | Path,
              keep_parent: bool = True) -> list[str]:
    """List all directories at a given root

    Args:
        root (str or Path): Path to directory whose folders need to be listed
        prefix (bool, optional): If true, prepends the path to each result, otherwise
            only returns the name of the directories found
    """
    root = Path(root).expanduser()
    directories = [str(p) for p in root.iterdir() if p.is_dir()]
    if not keep_parent:
        directories = [os.path.basename(d) for d in directories]
    return directories


def list_files(root: str | Path,
               prefix: str | tuple[str, ...] = "",
               suffix: str | tuple[str, ...] = "",
               keep_parent: bool = True) -> list[str]:
    """List all files at a given root

    Args:
        root (str or Path): Path to directory whose files need to be listed
        prefix (str or tuple): Prefix of the file names to match,
            Defaults to empty string '""'.
        suffix (str or tuple): Suffix of the files to match, e.g. ".png" or
            (".jpg", ".png").
            Defaults to empty string '""'.
        keep_parent (bool, optional): If true, prepends the parent path to each
            result, otherwise only returns the name of the files found.
            Defaults to False.
    """
    root = Path(root).expanduser()
    files = [
        str(p) for p in root.iterdir() if p.is_file()
        and p.name.startswith(prefix) and p.name.endswith(suffix)
    ]

    if not keep_parent:
        files = [os.path.basename(f) for f in files]

    return files


def _extract_tar(from_path: str | Path, to_path: str | Path,
                 compression: str | None) -> None:
    with tarfile.open(from_path,
                      f"r:{compression[1:]}" if compression else "r") as tar:
        tar.extractall(to_path)


_ZIP_COMPRESSION_MAP: dict[str, int] = {
    ".bz2": zipfile.ZIP_BZIP2,
    ".xz": zipfile.ZIP_LZMA,
}


def _extract_zip(from_path: str | Path, to_path: str | Path,
                 compression: str | None) -> None:
    with zipfile.ZipFile(from_path,
                         "r",
                         compression=_ZIP_COMPRESSION_MAP[compression]
                         if compression else zipfile.ZIP_STORED) as zip:
        zip.extractall(to_path)


_ARCHIVE_EXTRACTORS: dict[str, Callable[[str, str, str | None], None]] = {
    ".tar": _extract_tar,
    ".zip": _extract_zip,
}
_COMPRESSED_FILE_OPENERS: dict[str, Callable[..., IO]] = {
    ".bz2": bz2.open,
    ".gz": gzip.open,
    ".xz": lzma.open,
}
_FILE_TYPE_ALIASES: dict[str, tuple[str | None, str | None]] = {
    ".tbz": (".tar", ".bz2"),
    ".tbz2": (".tar", ".bz2"),
    ".tgz": (".tar", ".gz"),
}


def _detect_file_type(file: str | Path) -> tuple[str, str | None, str | None]:
    """Detect the archive type and/or compression of a file.

    Args:
        file (str, Path): the filename

    Returns:
        (tuple): tuple of suffix, archive type, and compression

    Raises:
        RuntimeError: if file has no suffix or suffix is not supported
    """
    suffixes = Path(file).suffixes
    if not suffixes:
        raise RuntimeError(
            f"File '{file}' has no suffixes that could be used to detect the archive type and compression."
        )
    suffix = suffixes[-1]

    # check if the suffix is a known alias
    if suffix in _FILE_TYPE_ALIASES:
        return (suffix, *_FILE_TYPE_ALIASES[suffix])

    # check if the suffix is an archive type
    if suffix in _ARCHIVE_EXTRACTORS:
        return suffix, suffix, None

    # check if the suffix is a compression
    if suffix in _COMPRESSED_FILE_OPENERS:
        # check for suffix hierarchy
        if len(suffixes) > 1:
            suffix2 = suffixes[-2]

            # check if the suffix2 is an archive type
            if suffix2 in _ARCHIVE_EXTRACTORS:
                return suffix2 + suffix, suffix2, suffix

        return suffix, None, suffix

    valid_suffixes = sorted(
        set(_FILE_TYPE_ALIASES)
        | set(_ARCHIVE_EXTRACTORS)
        | set(_COMPRESSED_FILE_OPENERS))
    raise RuntimeError(
        f"Unknown compression or archive type: '{suffix}'.\nKnown suffixes are: '{valid_suffixes}'."
    )


def _decompress(from_path: Path | str,
                to_path: Path | str | None = None,
                remove_finished: bool = False) -> str:
    r"""Decompress a file.

    The compression is automatically detected from the file name.

    Args:
        from_path (str or Path): Path to the file to be decompressed.
        to_path (str Path): Path to the decompressed file. If omitted, ``from_path`` without compression extension is used.
        remove_finished (bool): If ``True``, remove the file after the extraction.

    Returns:
        (str): Path to the decompressed file.
    """
    suffix, archive_type, compression = _detect_file_type(from_path)
    if not compression:
        raise RuntimeError(
            f"Couldn't detect a compression from suffix {suffix}.")

    if to_path is None:
        to_path = str(from_path).replace(
            suffix, archive_type if archive_type is not None else "")

    compressed_file_opener = _COMPRESSED_FILE_OPENERS[compression]

    with compressed_file_opener(from_path, "rb") as rfh, open(to_path,
                                                              "wb") as wfh:
        wfh.write(rfh.read())

    if remove_finished:
        os.remove(from_path)

    return str(to_path)


def extract_archive(from_path: str | Path,
                    to_path: str | Path | None = None,
                    remove_finished: bool = False) -> str:
    """Extract an archive.

    The archive type and a possible compression is automatically detected from
    the file name. If the file is compressed but not an archive the call is
    dispatched to :func:`decompress`.

    Args:
        from_path (str, Path): Path to the file to be extracted.
        to_path (str, Path): Path to the directory the file will be extracted to.
            If omitted, the directory of the archive file is used.
        remove_finished (bool): If ``True``, remove the file after the extraction.

    Returns:
        (str): Path to the directory the file was extracted to.
    """
    from_path = Path(from_path)

    if to_path is None:
        to_path = from_path.parent
    else:
        to_path = Path(to_path)

    suffix, archive_type, compression = _detect_file_type(from_path)
    if not archive_type:
        return _decompress(
            from_path,
            to_path / from_path.name.replace(suffix, ""),
            remove_finished=remove_finished,
        )

    extractor = _ARCHIVE_EXTRACTORS[archive_type]

    extractor(str(from_path), str(to_path), compression)
    if remove_finished:
        from_path.unlink()

    return str(to_path)

def extract_file_matching_pattern(archive: zipfile.ZipFile, prefix: str, suffix: str, extract_dir: Path, out_filename: str|None = None):
    """Extract a file matching a pattern from an archive and place it in the extraction directory under the given filename.

    Args:
        archive(zipfile.ZipFile): Archive from which to extract the file
        prefix(str): Prefix to match in the filename. Pass empty string for no prefix.
        suffix(str): Suffix to match in the filename. Pass empty string for no suffix.
        extract_dir(Path): Path to the folder where to store the extracted file
        out_filename(str): Name to assign to the extracted file.

    Examples:
        >>> extract_file_matching_pattern(zipfile.ZipFile("archive.zip"), "", ".txt", ".", "results.txt")
        """
    files: list[str] = [x.filename for x in archive.filelist]
    file_to_extract = list(filter(lambda x: x.startswith(prefix) and x.endswith(suffix), files)).pop()
    archive.extract(file_to_extract, extract_dir)
    if out_filename is not None:
        os.rename(extract_dir / file_to_extract, extract_dir / out_filename)


def download_and_extract_archive(
    url: str,
    download_root: str | Path,
    extract_root: str | Path | None = None,
    filename: str | None = None,
    md5: str | None = None,
    remove_finished: bool = False,
) -> None:
    """Download a file from url and extract it

       This method is a wrapper of `download_url` and `extract_archive` methods.

    Args:
        url (str): URL to download file from
        download_root (str or Path): Path to the directory to place downloaded
            file in
        extract_root (str or Path, optional): Path to the directory the file
            will be extracted to. If omitted, the `download_root` is used.
        filename (str, optional): Name to save the downloaded file under.
            If None, use the basename of the URL
        md5 (str, optional): MD5 checksum of the download. If None, do not check
        remove_finished (bool, optional): If `True`, remove the downloaded file
             after the extraction. Defaults to False.
    """

    download_root = Path(download_root).expanduser()
    if extract_root is None:
        extract_root = download_root
    else:
        extract_root = Path(extract_root)
    if not filename:
        filename = Path(url).name

    download_url(url, download_root, filename, md5)

    archive = download_root / filename
    print(f"Extracting {archive} to {extract_root}")
    extract_archive(archive, extract_root, remove_finished)
