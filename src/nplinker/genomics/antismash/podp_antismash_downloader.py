from __future__ import annotations
import json
import logging
import time
import warnings
from collections.abc import Mapping
from collections.abc import Sequence
from os import PathLike
from pathlib import Path
from jsonschema import validate
from nplinker.defaults import GENOME_STATUS_FILENAME
from nplinker.genomics.antismash import antismash_job_is_done
from nplinker.genomics.antismash import download_and_extract_from_antismash_api
from nplinker.genomics.antismash import download_and_extract_from_antismash_db
from nplinker.genomics.antismash import download_and_extract_ncbi_genome
from nplinker.genomics.antismash import extract_antismash_data
from nplinker.genomics.antismash import resolve_genome_accession
from nplinker.genomics.antismash import submit_antismash_job
from nplinker.schemas import GENOME_STATUS_SCHEMA


logger = logging.getLogger(__name__)

JGI_GENOME_LOOKUP_URL = (
    "https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid={}"
)
USER_AGENT = "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:86.0) Gecko/20100101 Firefox/86.0"


class GenomeStatus:
    """Class to represent the status of a single genome.

    The status of genomes is tracked in the file
    [GENOME_STATUS_FILENAME][nplinker.defaults.GENOME_STATUS_FILENAME].
    """

    def __init__(
        self,
        original_id: str,
        resolved_id: str = "",
        failed_previously: bool = False,
        bgc_path: str = "",
    ):
        """Initialize a GenomeStatus object for the given genome.

        Args:
            original_id: The original ID of the genome.
            resolved_id: The resolved genome ID of the genome. Defaults to "".
            failed_previously: Indicates whether a previous attempt to get BGC data
                for the genome has failed. Defaults to False.
            bgc_path: The path to the downloaded BGC file for
                the genome. Defaults to "".
        """
        self.original_id = original_id
        self.resolved_id = "" if resolved_id == "None" else resolved_id
        self.failed_previously = failed_previously
        self.bgc_path = bgc_path

    @staticmethod
    def read_json(file: str | PathLike) -> dict[str, "GenomeStatus"]:
        """Get a dict of GenomeStatus objects by loading given genome status file.

        Note that an empty dict is returned if the given file doesn't exist.

        Args:
            file: Path to genome status file.

        Returns:
            Dict keys are genome original id and values are GenomeStatus
                objects. An empty dict is returned if the given file doesn't exist.
        """
        genome_status_dict = {}
        if Path(file).exists():
            with open(file, "r") as f:
                data = json.load(f)

            # validate json data before using it
            validate(data, schema=GENOME_STATUS_SCHEMA)

            genome_status_dict = {
                gs["original_id"]: GenomeStatus(**gs) for gs in data["genome_status"]
            }
        return genome_status_dict

    @staticmethod
    def to_json(
        genome_status_dict: Mapping[str, "GenomeStatus"], file: str | PathLike | None = None
    ) -> str | None:
        """Convert the genome status dictionary to a JSON string.

        If a file path is provided, the JSON string is written to the file. If
        the file already exists, it is overwritten.

        Args:
            genome_status_dict: A dictionary of genome
                status objects. The keys are the original genome IDs and the values
                are GenomeStatus objects.
            file: The path to the output JSON file.
                If None, the JSON string is returned but not written to a file.

        Returns:
            The JSON string if `file` is None, otherwise None.
        """
        gs_list = [gs._to_dict() for gs in genome_status_dict.values()]
        json_data = {"genome_status": gs_list, "version": "1.0"}

        # validate json object before dumping
        validate(json_data, schema=GENOME_STATUS_SCHEMA)

        if file is not None:
            with open(file, "w") as f:
                json.dump(json_data, f)
            return None
        return json.dumps(json_data)

    def _to_dict(self) -> dict:
        """Convert the GenomeStatus object to a dict."""
        return {
            "original_id": self.original_id,
            "resolved_id": self.resolved_id,
            "failed_previously": self.failed_previously,
            "bgc_path": self.bgc_path,
        }


def podp_download_and_extract_antismash_data(
    genome_records: Sequence[Mapping[str, Mapping[str, str]]],
    project_download_root: str | PathLike,
    project_extract_root: str | PathLike,
):
    """Download and extract antiSMASH BGC archive for the given genome records.

    Args:
        genome_records: list of dicts representing genome records.

            The dict of each genome record contains a key of genome ID with a value
            of another dict containing information about genome type, label and
            accession ids (RefSeq, GenBank, and/or JGI).
        project_download_root: Path to the directory to place
            downloaded archive in.
        project_extract_root: Path to the directory downloaded archive will be extracted to.

            Note that an `antismash` directory will be created in the specified
            `extract_root` if it doesn't exist. The files will be extracted to
            `<extract_root>/antismash/<antismash_id>` directory.

    Warnings:
        UserWarning: when no antiSMASH data is found for some genomes.
    """
    if not Path(project_download_root).exists():
        # otherwise in case of failed first download, the folder doesn't exist and
        # genome_status_file can't be written
        Path(project_download_root).mkdir(parents=True, exist_ok=True)

    gs_file = Path(project_download_root, GENOME_STATUS_FILENAME)
    gs_dict = GenomeStatus.read_json(gs_file)

    for i, genome_record in enumerate(genome_records):
        logger.info(
            f"Getting antismash BGC data for genome record {i + 1} of {len(genome_records)}."
        )

        # get the best available genome ID from the dict
        original_genome_id = get_best_available_genome_id(genome_record["genome_ID"])
        if not original_genome_id:
            logger.warning(f"Skipping invalid genome record: {genome_record}")
            continue
        # Retrieve or initialize the GenomeStatus object for the genome ID
        gs = gs_dict.setdefault(original_genome_id, GenomeStatus(original_genome_id))

        # Check if genomes already have antiSMASH BGC data
        if gs.bgc_path and Path(gs.bgc_path).exists():
            logger.info(
                f"antiSMASH BGC data for genome ID {original_genome_id} already downloaded to "
                f"{gs.bgc_path}"
            )
            try:
                process_existing_antismash_data(gs, project_extract_root)
                continue
            except Exception as e:
                logger.warning(
                    "Failed to process existing antiSMASH BGC data for genome ID "
                    f"{original_genome_id}. Error: {e}"
                )
        gs.bgc_path = ""  # Reset bgc path

        # Check if a previous attempt to get bgc data has failed
        if gs.failed_previously:
            logger.info(f"Genome ID {original_genome_id} skipped due to previous failed attempt")
            continue

        # resolve genome ID
        try:
            gs.resolved_id = resolve_genome_accession(genome_record["genome_ID"])
        except Exception as e:
            logger.warning(f"Failed to resolve genome ID {gs.original_id}. Error: {e}")
            gs.failed_previously = True
            continue

        # retrieve antismash BGC data from antiSMASH-DB
        try:
            retrieve_antismash_db_data(gs, project_download_root, project_extract_root)
            logger.info(
                f"antiSMASH BGC data for genome ID {gs.original_id} is downloaded and extracted"
            )
            continue
        except Exception as e:
            logger.info(
                f"Unable to retrieve BGC data from antiSMASH-DB for genome ID {gs.original_id}. "
                f"Error: {e}"
            )

        # retrieve antismash BGC by submitting antismash job via API
        try:
            logger.info(
                "Downloading genome assembly from NCBI and submitting antiSMASH job for "
                f"genome ID {gs.original_id}."
            )
            genome_path = download_and_extract_ncbi_genome(
                gs.resolved_id, project_download_root, project_extract_root
            )
            job_id = submit_antismash_job(genome_path)
            logger.info(f"Waiting for antiSMASH job {job_id} to complete.")
            while antismash_job_is_done(job_id) is False:
                time.sleep(15)
            retrieve_antismash_job_data(job_id, gs, project_download_root, project_extract_root)
            logger.info(
                f"antiSMASH BGC data for genome ID {gs.original_id} is downloaded and extracted"
            )
            continue
        except Exception as e:
            logger.info(
                f"Unable to retrieve BGC data via antiSMASH API for genome ID {gs.original_id}. "
                f"Error: {e}"
            )

        if gs.bgc_path == "":
            logger.warning(f"Failed to retrieve BGC data for genome ID {gs.original_id}.")
            gs.failed_previously = True

    # raise and log warning for failed downloads
    failed_ids = [gs.original_id for gs in gs_dict.values() if not gs.bgc_path]
    if failed_ids:
        warning_message = (
            f"Failed to download antiSMASH data for the following genome IDs: {failed_ids}"
        )
        logger.warning(warning_message)
        warnings.warn(warning_message, UserWarning)

    # save updated genome status to json file
    GenomeStatus.to_json(gs_dict, gs_file)

    if len(failed_ids) == len(genome_records):
        raise ValueError("No antiSMASH data found for any genome")


def get_best_available_genome_id(genome_id_data: Mapping[str, str]) -> str | None:
    """Get the best available ID from genome_id_data dict.

    Args:
        genome_id_data: dictionary containing information for each genome record present.

    Returns:
        ID for the genome, if present, otherwise None.
    """
    if "RefSeq_accession" in genome_id_data:
        best_id = genome_id_data["RefSeq_accession"]
    elif "GenBank_accession" in genome_id_data:
        best_id = genome_id_data["GenBank_accession"]
    elif "JGI_Genome_ID" in genome_id_data:
        best_id = genome_id_data["JGI_Genome_ID"]
    else:
        best_id = None

    if best_id is None or len(best_id) == 0:
        logger.warning(f"Failed to get valid genome ID in genome data: {genome_id_data}")
        return None
    return best_id


def process_existing_antismash_data(gs_obj: GenomeStatus, extract_root: str | PathLike) -> None:
    """Processes already downloaded antiSMASH BGC data archive.

    This function ensures that the antiSMASH data archive associated with a given genomic sequence
    object is properly extracted into a specified directory. If the data has already been extracted,
    the function skips the extraction process.

    Args:
        gs_obj: An object representing a genomic sequence, which contains the path
                to the antiSMASH BGC data (accessible via `gs_obj.bgc_path`) and
                an original identifier (`gs_obj.original_id`).
        extract_root: The root directory where the antiSMASH data should be extracted.

    Raises:
        Any exceptions raised by the `extract_antismash_data` function if the extraction fails.
    """
    antismash_id = Path(gs_obj.bgc_path).stem
    extract_path = Path(extract_root, "antismash", antismash_id)
    completed_marker = extract_path / "completed"

    # Check if archive is already successfully extracted
    if completed_marker.exists():
        logger.info(
            f"antiSMASH BGC data for {gs_obj.original_id} already extracted at {extract_path}."
        )
        return

    extract_antismash_data(gs_obj.bgc_path, extract_root, antismash_id)
    completed_marker.touch(exist_ok=True)


def retrieve_antismash_db_data(
    genome_status: GenomeStatus, download_root: str | PathLike, extract_root: str | PathLike
) -> None:
    """Retrieve antiSMASH database data for a given genome and update its status.

    This function downloads and extracts antiSMASH data for a genome identified
    by its resolved genome ID. It updates the `genome_status` object with the
    path to the downloaded data or sets it to an empty string if an error occurs.

    Args:
        genome_status (GenomeStatus): An object representing the genome's status,
            including its resolved genome ID and BGC path.
        download_root (str | PathLike): The root directory where the antiSMASH
            data will be downloaded.
        extract_root (str | PathLike): The root directory where the antiSMASH
            data will be extracted.

    Raises:
        Exception: If an error occurs during the download or extraction process.
    """
    if not genome_status.resolved_id.startswith("GCF_"):
        raise ValueError(
            f"Resolved genome ID '{genome_status.resolved_id}' is not a valid RefSeq assembly and "
            "antiSMASH-DB only contains results for RefSeq assemblies."
        )

    antismash_id = genome_status.resolved_id
    extract_path = Path(extract_root, "antismash", antismash_id)
    download_path = Path(download_root, f"{antismash_id}.zip").absolute()

    download_and_extract_from_antismash_db(antismash_id, download_root, extract_root)
    Path.touch(extract_path / "completed", exist_ok=True)
    genome_status.bgc_path = str(download_path)


def retrieve_antismash_job_data(
    job_id: str,
    genome_status: GenomeStatus,
    download_root: str | PathLike,
    extract_root: str | PathLike,
) -> None:
    """Retrieve antiSMASH API data for a given genome and update its status.

    This function downloads and extracts antiSMASH data for a genome identified
    by its resolved genome ID. It updates the `genome_status` object with the
    path to the downloaded data or sets it to an empty string if an error occurs.

    Args:
        job_id (str): The job ID for the antiSMASH API job.
        genome_status (GenomeStatus): An object representing the genome's status,
            including its resolved genome ID and BGC path.
        download_root (str | PathLike): The root directory where the antiSMASH
            data will be downloaded.
        extract_root (str | PathLike): The root directory where the antiSMASH
            data will be extracted.

    Raises:
        Exception: If an error occurs during the download or extraction process.
    """
    antismash_id = genome_status.resolved_id
    extract_path = Path(extract_root, "antismash", antismash_id)
    download_path = Path(download_root, f"{antismash_id}.zip").absolute()

    download_and_extract_from_antismash_api(job_id, antismash_id, download_root, extract_root)
    Path.touch(extract_path / "completed", exist_ok=True)
    genome_status.bgc_path = str(download_path)
