from enum import Enum
import os
from os import PathLike
from pathlib import Path
import zipfile
from bs4 import BeautifulSoup
from bs4 import Tag
import requests
from nplinker.utils import get_headers


# TODO: add description
class GNPSFormat(Enum):
    Unknown = 0
    AllFiles = 1
    UniqueFiles = 2
    FBMN = 3

GNPS_TASK_URL = 'https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task={}'


def gnps_format_from_file_mapping(file: str | PathLike) -> GNPSFormat:
    """Detect GNPS format from the given file mapping file.

    Different GNPS workflows has different GNSP file mapping file:
    - METABOLOMICS-SNETS workflow: the .tsv file under folder "clusterinfosummarygroup_attributes_withIDs_withcomponentID"
    - METABOLOMICS-SNETS-V2 workflow: the .clustersummary file (tsv) under folder "clusterinfosummarygroup_attributes_withIDs_withcomponentID"
    - FEATURE-BASED-MOLECULAR-NETWORKING workflow: the .csv file under folder "quantification_table"

    Args:
        file(str | PathLike): Path to the file to peek the format for.

    Returns:
        GNPSFormat: GNPS format identified in the file.
    """
    headers = get_headers(file)
    if 'AllFiles' in headers:
        return GNPSFormat.AllFiles
    if 'UniqueFileSources' in headers:
        return GNPSFormat.UniqueFiles
    if 'row ID' in headers:
        return GNPSFormat.FBMN
    return GNPSFormat.Unknown


def gnps_format_from_task_id(task_id: str) -> GNPSFormat:
    """Detect GNPS format for the given task id.

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
    workflow_tag: Tag = list(filter(lambda x: x.contents == ['Workflow'],
                                    tags))[0]
    workflow_format_tag: Tag = workflow_tag.parent.contents[3]
    workflow_format = workflow_format_tag.contents[0].strip()

    if workflow_format == "FEATURE-BASED-MOLECULAR-NETWORKING":
        return GNPSFormat.FBMN
    if workflow_format == "METABOLOMICS-SNETS-V2":
        return GNPSFormat.UniqueFiles
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
