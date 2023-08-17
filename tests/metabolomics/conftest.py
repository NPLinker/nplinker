from os import PathLike
import shutil
import httpx
import pytest
from nplinker.metabolomics.gnps.gnps_format import GNPSFormat
from nplinker.utils import extract_archive
from .. import GNPS_DATA_DIR


@pytest.fixture(scope="session", autouse=True)
def gnps_website_is_down():
    gnps_url = "https://gnps.ucsd.edu"
    try:
        _ = httpx.get(gnps_url)
        return False
    except httpx.HTTPError:
        return True


@pytest.fixture(scope="session", autouse=True)
def gnps_zip_files() -> dict[GNPSFormat, PathLike]:
    """Get the GNPS zip archives as a dictionary.

    The dict keys are the workflow short names taken from the GNPSFormat enum.
    The dict values are the paths to the zip archives.

    The GNPS zip archives are not complete and contains only the files required
    for NPLinker. You can download the archives from the following links:
    - https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=c22f44b14a3d450eb836d607cb9521bb
    - https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=189e8bf16af145758b0a900f1c44ff4a
    - https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=92036537c21b44c29e509291e53f6382
    """
    return {
        GNPSFormat.SNETS: GNPS_DATA_DIR /
        "ProteoSAFe-METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra.zip",
        GNPSFormat.SNETSV2: GNPS_DATA_DIR /
        "ProteoSAFe-METABOLOMICS-SNETS-V2-189e8bf1-download_clustered_spectra.zip",
        GNPSFormat.FBMN: GNPS_DATA_DIR /
        "ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-92036537-download_cytoscape_data.zip",
        GNPSFormat.Unknown: GNPS_DATA_DIR / "ProteoSAFe-Unknown.zip"
    }


@pytest.fixture(scope="session", autouse=True)
def prepare_data(gnps_zip_files):
    """Extract GNPS zip archive to a temporary directory and delete it after the test run.

    The temporary directory is named after the workflow, e.g. "SNETS", "SNETSV2", "FBMN".
    Note that these directory names are also used in other fixtures.
    """
    for workflow, zip_file in gnps_zip_files.items():
        extract_archive(zip_file, GNPS_DATA_DIR / workflow.name)
    yield
    for workflow in gnps_zip_files:
        shutil.rmtree(GNPS_DATA_DIR / workflow.name)


@pytest.fixture(scope="session")
def gnps_file_mappings_files() -> dict[GNPSFormat, PathLike]:
    return {
        GNPSFormat.SNETS:
        GNPS_DATA_DIR / GNPSFormat.SNETS.name /
        "clusterinfosummarygroup_attributes_withIDs_withcomponentID" /
        "d69356c8e5044c2a9fef3dd2a2f991e1.tsv",
        GNPSFormat.SNETSV2:
        GNPS_DATA_DIR / GNPSFormat.SNETSV2.name /
        "clusterinfosummarygroup_attributes_withIDs_withcomponentID" /
        "16f782af01bc4f50a23ed163566072f9.clustersummary",
        GNPSFormat.FBMN:
        GNPS_DATA_DIR / GNPSFormat.FBMN.name /
        "quantification_table_reformatted" /
        "1a12f6fbd2ca4e099ec56bdaea56368f.csv"
    }


@pytest.fixture(scope="session")
def gnps_spectra_files() -> dict[GNPSFormat, PathLike]:
    return {
        GNPSFormat.SNETS:
        GNPS_DATA_DIR / GNPSFormat.SNETS.name /
        "METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra-main.mgf",
        GNPSFormat.SNETSV2:
        GNPS_DATA_DIR / GNPSFormat.SNETSV2.name /
        "METABOLOMICS-SNETS-V2-189e8bf1-download_clustered_spectra-main.mgf",
        GNPSFormat.FBMN:
        GNPS_DATA_DIR / GNPSFormat.FBMN.name / "spectra" / "specs_ms.mgf"
    }


@pytest.fixture(scope="session")
def gnps_mf_files() -> dict[GNPSFormat, PathLike]:
    return {
        GNPSFormat.SNETS:
        GNPS_DATA_DIR / GNPSFormat.SNETS.name / "networkedges_selfloop" /
        "6da5be36f5b14e878860167fa07004d6.pairsinfo",
        GNPSFormat.SNETSV2:
        GNPS_DATA_DIR / GNPSFormat.SNETSV2.name / "networkedges_selfloop" /
        "06dd31e28bb547ba852859219db9298c..selfloop",
        GNPSFormat.FBMN:
        GNPS_DATA_DIR / GNPSFormat.FBMN.name / "networkedges_selfloop" /
        "c74fec018736475483e9c8b05e230cce..selfloop"
    }


@pytest.fixture(scope="session")
def gnps_annotations_files() -> dict[GNPSFormat, PathLike]:
    return {
        GNPSFormat.SNETS:
        GNPS_DATA_DIR / GNPSFormat.SNETS.name / "result_specnets_DB" /
        "885e4c5485ba42569e4876d1fe90d759.tsv",
        GNPSFormat.SNETSV2:
        GNPS_DATA_DIR / GNPSFormat.SNETSV2.name / "result_specnets_DB" /
        "017fadadf6744c10b5d39f109e1438dc.tsv",
        GNPSFormat.FBMN:
        GNPS_DATA_DIR / GNPSFormat.FBMN.name / "DB_result" /
        "7dc5b46b50d94246a1de12ef485d0f75.tsv"
    }
