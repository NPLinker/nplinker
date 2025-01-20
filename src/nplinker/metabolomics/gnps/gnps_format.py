from __future__ import annotations
import tarfile
import zipfile
from enum import Enum
from enum import unique
from os import PathLike
from pathlib import Path
import httpx
import yaml
from bs4 import BeautifulSoup


@unique
class GNPSFormat(Enum):
    """Enum class for GNPS formats or workflows.

    ??? info "Concept"
        [GNPS data][gnps-data]

    The name of the enum is a short name for the workflow, and the value of the enum is the workflow
    name used on the GNPS website.
    """

    # Format: ShortName = "GNPSWorkflowName"
    # For GNPS1
    SNETS = "METABOLOMICS-SNETS"
    SNETSV2 = "METABOLOMICS-SNETS-V2"
    FBMN = "FEATURE-BASED-MOLECULAR-NETWORKING"
    # For GNPS2
    GNPS2CN = "classical_networking_workflow"
    GNPS2FBMN = "feature_based_molecular_networking_workflow"
    # Unknown format
    Unknown = "Unknown-GNPS-Workflow"


def gnps_format_from_gnps1_task_id(task_id: str) -> GNPSFormat:
    """Detect GNPS format or workflow for the given GNPS1 task id.

    GNPS1 tasks are those generated on the platform https://gnps.ucsd.edu.

    Args:
        task_id: GNPS1 task id.

    Returns:
        The format identified in the task.

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
    gnps1_task_url = "https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task={}"
    task_html = httpx.get(gnps1_task_url.format(task_id))
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


def gnps_format_from_archive(file: str | PathLike) -> GNPSFormat:
    """Detect GNPS format or workflow from GNPS archive file.

    GNPS archive files can be in two formats: GNPS1 (.zip) and GNPS2 (.tar).

    For GNPS1 data, the detection of workflow format is based on the filename of the zip archive and
    the names of the files contained in the zip archive.

    For GNPS2 data, the workflow format is taken from the `submission_parameters.yaml` file in the
    tar archive, which has a key `workflowname`.

    Args:
        file: Path to the GNPS archive file.

    Returns:
        The format identified in the GNPS archive file.

    Examples:
        >>> gnps_format_from_archive("ProteoSAFe-METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra.zip")
        <GNPSFormat.SNETS: 'METABOLOMICS-SNETS'>
        >>> gnps_format_from_archive("ProteoSAFe-METABOLOMICS-SNETS-V2-189e8bf1-download_clustered_spectra.zip")
        <GNPSFormat.SNETSV2: 'METABOLOMICS-SNETS-V2'>
        >>> gnps_format_from_archive("ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-672d0a53-download_cytoscape_data.zip")
        <GNPSFormat.FBMN: 'FEATURE-BASED-MOLECULAR-NETWORKING'>
        >>> gnps_format_from_archive("206a7b40b7ed41c1ae6b4fbd2def3636.tar")
        <GNPSFormat.GNPS2CN: 'classical_networking_workflow'>
        >>> gnps_format_from_archive("2014f321d72542afb5216c932e0d5079.tar")
        <GNPSFormat.GNPS2FBMN: 'feature_based_molecular_networking_workflow'>
    """
    file = Path(file)
    suffix = file.suffix
    if suffix == ".zip":
        return _gnps_format_from_archive_gnps1(file)
    if suffix == ".tar":
        return _gnps_format_from_archive_gnps2(file)
    return GNPSFormat.Unknown


def _gnps_format_from_archive_gnps1(file: Path) -> GNPSFormat:
    """Detect GNPS format from GNPS1 archive file."""
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


def _gnps_format_from_archive_gnps2(file: Path) -> GNPSFormat:
    """Detect GNPS format from GNPS2 archive file."""
    with tarfile.open(file, "r") as tar:
        try:
            submission_file = tar.extractfile("submission_parameters.yaml")
            if submission_file is None:
                return GNPSFormat.Unknown
            submission_params = yaml.safe_load(submission_file)
        except (KeyError, yaml.YAMLError):
            return GNPSFormat.Unknown

    workflow = submission_params.get("workflowname")

    if workflow == GNPSFormat.GNPS2FBMN.value:
        return GNPSFormat.GNPS2FBMN
    if workflow == GNPSFormat.GNPS2CN.value:
        return GNPSFormat.GNPS2CN
    return GNPSFormat.Unknown
