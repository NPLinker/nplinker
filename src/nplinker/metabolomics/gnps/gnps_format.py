from enum import Enum
import os
from os import PathLike
from pathlib import Path
import zipfile
from bs4 import BeautifulSoup
from bs4 import Tag
import requests
from nplinker.utils import get_headers


class GNPSFormat(Enum):
    Unknown = 0
    AllFiles = 1
    UniqueFiles = 2
    FBMN = 3

GNPS_TASK_URL = 'https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task={}'


def gnps_format_from_file_mapping(file: str | PathLike, has_quant_table: bool) -> GNPSFormat:
    """Peek GNPS file format for given file.

    TODO: #89 This should be rewritten to actually return the format always based on only the file and not include the quant table in it.

     Args:
        file(str | PathLike): Path to the file to peek the format for.
        has_quant_table (bool): If a quant table is present, do return GNPS_FORMAT_NEW_FBMN.

    Returns:
        GNPSFormat: GNPS format identified in the file.
    """

    headers: list[str] = get_headers(os.fspath(file))

    if headers is None:
        return GNPSFormat.Unknown

    # first, check for AllFiles
    if 'AllFiles' in headers:
        # this should be an old-style dataset like Crusemann, with a single .tsv file
        # containing all the necessary info. The AllFiles column should contain pairs
        # of mzXML filenames and scan numbers in this format:
        #   filename1:scannumber1###filename2:scannumber2###...
        return GNPSFormat.AllFiles
    elif 'UniqueFileSources' in headers:
        # this is a slightly newer-style dataset, e.g. MSV000084771 on the platform
        # it still just has a single .tsv file, but AllFiles is apparently replaced
        # with a UniqueFileSources column. There is also a UniqueFileSourcesCount
        # column which just indicates the number of entries in the UniqueFileSources
        # column. If there are multiple entries the delimiter is a | character
        return GNPSFormat.UniqueFiles
    elif has_quant_table:
        # if there is no AllFiles/UniqueFileSources, but we DO have a quantification
        # table file, that should indicate a new-style dataset like Carnegie
        # TODO check for the required header columns here too
        return GNPSFormat.FBMN
    elif len(list(filter(lambda x: "Peak area" in x, headers))) > 1:
        return GNPSFormat.FBMN
    else:
        # if we don't match any of the above cases then it's not a recognised format
        return GNPSFormat.Unknown


def gnps_format_from_task_id(task_id: str) -> GNPSFormat:
    """Detect the GNPS format for the given task id.

    The http request has a timeout of 5 seconds. If the request fails,
    an ReadTimeout exception is raised. This is to prevent the program
    from hanging indefinitely when the GNPS server is down.

    Args:
        task_id(str): GNPS task id.

    Returns:
        GNPSFormat: the format used in the GNPS workflow invocation.

    Examples:
        >>> gnps_format_from_task_id("92036537c21b44c29e509291e53f6382")
    """
    task_html = requests.get(GNPS_TASK_URL.format(task_id), timeout=5)
    soup = BeautifulSoup(task_html.text, features="html.parser")
    tags = soup.find_all('th')
    workflow_tag: Tag = list(filter(lambda x: x.contents == ['Workflow'], tags))[0]
    workflow_format_tag: Tag = workflow_tag.parent.contents[3]
    workflow_format = workflow_format_tag.contents[0].strip()

    if workflow_format == "FEATURE-BASED-MOLECULAR-NETWORKING":
        return GNPSFormat.FBMN
    if workflow_format == "METABOLOMICS-SNETS":
        return GNPSFormat.AllFiles
    return GNPSFormat.Unknown


def gnps_format_from_archive(zip_file: str | PathLike) -> GNPSFormat:
    """Detect GNPS format from a downloaded zip archive.

    The detection is based on the filename of the zip file and the names of the
    files contained in the zip file.

    Args:
        zip_file: Path to the downloaded zip file.

    Returns:
        GNPSFormat: the format used in the GNPS workflow invocation.

    Examples:
        >>> gnps_format_from_archive("tests/data/ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-92036537-download_cytoscape_data.zip")
     """
    file = Path(zip_file)
    # Guess the format from the filename of the zip file
    if "FEATURE-BASED-MOLECULAR-NETWORKING" in file.name:
        return GNPSFormat.FBMN
    if "METABOLOMICS-SNETS-V2" in file.name:
        return GNPSFormat.UniqueFiles
    if "METABOLOMICS-SNETS" in file.name:
        return GNPSFormat.AllFiles

    # Guess the format from the names of the files in the zip file
    with zipfile.ZipFile(file) as archive:
        filenames = archive.namelist()
    if any("FEATURE-BASED-MOLECULAR-NETWORKING" in x for x in filenames):
        return GNPSFormat.FBMN
    if any("METABOLOMICS-SNETS-V2" in x for x in filenames):
        return GNPSFormat.UniqueFiles
    if any("METABOLOMICS-SNETS" in x for x in filenames):
        return GNPSFormat.AllFiles

    return GNPSFormat.Unknown
