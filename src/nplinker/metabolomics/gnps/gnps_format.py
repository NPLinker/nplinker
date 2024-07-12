from __future__ import annotations
import re
import zipfile
from enum import Enum
from enum import unique
from os import PathLike
from pathlib import Path
import httpx
from bs4 import BeautifulSoup


GNPS_TASK_URL = "https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task={}"


@unique
class GNPSFormat(Enum):
    """Enum class for GNPS formats or workflows.

    ??? info "Concept"
        [GNPS data][gnps-data]

    The name of the enum is a short name for the workflow, and the value of the enum is the workflow
    name used on the GNPS website.
    """

    # Format: ShortName = "GNPSWorkflowName"
    SNETS = "METABOLOMICS-SNETS"
    SNETSV2 = "METABOLOMICS-SNETS-V2"
    FBMN = "FEATURE-BASED-MOLECULAR-NETWORKING"
    Unknown = "Unknown-GNPS-Workflow"


def gnps_format_from_task_id(task_id: str) -> GNPSFormat:
    """Detect GNPS format for the given task id.

    Args:
        task_id: GNPS task id.

    Returns:
        The format identified in the GNPS task.

    Examples:
        >>> gnps_format_from_task_id("c22f44b14a3d450eb836d607cb9521bb")
        <GNPSFormat.SNETS: 'METABOLOMICS-SNETS'>
        >>> gnps_format_from_task_id("189e8bf16af145758b0a900f1c44ff4a")
        <GNPSFormat.SNETSV2: 'METABOLOMICS-SNETS-V2'>
        >>> gnps_format_from_task_id("92036537c21b44c29e509291e53f6382")
        <GNPSFormat.FBMN: 'FEATURE-BASED-MOLECULAR-NETWORKING'>
        >>> gnps_format_from_task_id("0ad6535e34d449788f297e712f43068a")
        <GNPSFormat.Unknown: 'Unknown-GNPS-Workflow'>
    """
    task_html = httpx.get(GNPS_TASK_URL.format(task_id))
    soup = BeautifulSoup(task_html.text, features="html.parser")
    try:
        # find the td tag that follows the th tag containing 'Workflow'
        workflow_tag = soup.find("th", string="Workflow").find_next_sibling("td")  # type: ignore
        workflow_format = workflow_tag.contents[0].strip()  # type: ignore
    except AttributeError:
        return GNPSFormat.Unknown

    if workflow_format == GNPSFormat.FBMN.value:
        return GNPSFormat.FBMN
    if workflow_format == GNPSFormat.SNETSV2.value:
        return GNPSFormat.SNETSV2
    if workflow_format == GNPSFormat.SNETS.value:
        return GNPSFormat.SNETS
    return GNPSFormat.Unknown


def gnps_format_from_archive(zip_file: str | PathLike) -> GNPSFormat:
    """Detect GNPS format from GNPS zip archive.

    The detection is based on the filename of the zip file and the names of the
    files contained in the zip file.

    Args:
        zip_file: Path to the GNPS zip file.

    Returns:
        The format identified in the GNPS zip file.

    Examples:
        >>> gnps_format_from_archive("ProteoSAFe-METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra.zip")
        <GNPSFormat.SNETS: 'METABOLOMICS-SNETS'>
        >>> gnps_format_from_archive("ProteoSAFe-METABOLOMICS-SNETS-V2-189e8bf1-download_clustered_spectra.zip")
        <GNPSFormat.SNETSV2: 'METABOLOMICS-SNETS-V2'>
        >>> gnps_format_from_archive("ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-672d0a53-download_cytoscape_data.zip")
        <GNPSFormat.FBMN: 'FEATURE-BASED-MOLECULAR-NETWORKING'>
    """
    file = Path(zip_file)
    # Guess the format from the filename of the zip file
    if GNPSFormat.FBMN.value in file.name:
        return GNPSFormat.FBMN
    # the order of the if statements matters for the following two
    if GNPSFormat.SNETSV2.value in file.name:
        return GNPSFormat.SNETSV2
    if GNPSFormat.SNETS.value in file.name:
        return GNPSFormat.SNETS

    # Guess the format from the names of the files in the zip file
    with zipfile.ZipFile(file) as archive:
        filenames = archive.namelist()
    if any(GNPSFormat.FBMN.value in x for x in filenames):
        return GNPSFormat.FBMN
    # the order of the if statements matters for the following two
    if any(GNPSFormat.SNETSV2.value in x for x in filenames):
        return GNPSFormat.SNETSV2
    if any(GNPSFormat.SNETS.value in x for x in filenames):
        return GNPSFormat.SNETS

    return GNPSFormat.Unknown


def gnps_format_from_file_mapping(file: str | PathLike) -> GNPSFormat:
    """Detect GNPS format from the given file mapping file.

    The GNPS file mapping file is located in different folders depending on the
    GNPS workflow. Here are the locations in corresponding GNPS zip archives:

    - `METABOLOMICS-SNETS` workflow: the `.tsv` file in the folder
        `clusterinfosummarygroup_attributes_withIDs_withcomponentID`
    - `METABOLOMICS-SNETS-V2` workflow: the `.clustersummary` file (tsv) in the folder
        `clusterinfosummarygroup_attributes_withIDs_withcomponentID`
    - `FEATURE-BASED-MOLECULAR-NETWORKING` workflow: the `.csv` file in the folder
        `quantification_table`

    Args:
        file: Path to the file to peek the format for.

    Returns:
        GNPS format identified in the file.
    """
    with open(file, "r") as f:
        header = f.readline().strip()

    if re.search(r"\bAllFiles\b", header):
        return GNPSFormat.SNETS
    if re.search(r"\bUniqueFileSources\b", header):
        return GNPSFormat.SNETSV2
    if re.search(r"\b{}\b".format(re.escape("row ID")), header):
        return GNPSFormat.FBMN
    return GNPSFormat.Unknown
