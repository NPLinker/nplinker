from __future__ import annotations
import json
import logging
import re
import time
import warnings
from collections.abc import Mapping
from collections.abc import Sequence
from os import PathLike
from pathlib import Path
import httpx
from bs4 import BeautifulSoup
from bs4 import NavigableString
from bs4 import Tag
from jsonschema import validate
from nplinker.defaults import GENOME_STATUS_FILENAME
from nplinker.genomics.antismash import download_and_extract_antismash_data
from nplinker.schemas import GENOME_STATUS_SCHEMA


logger = logging.getLogger(__name__)

NCBI_LOOKUP_URL = "https://www.ncbi.nlm.nih.gov/assembly/?term={}"
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
        resolved_refseq_id: str = "",
        resolve_attempted: bool = False,
        bgc_path: str = "",
    ):
        """Initialize a GenomeStatus object for the given genome.

        Args:
            original_id: The original ID of the genome.
            resolved_refseq_id: The resolved RefSeq ID of the
                genome. Defaults to "".
            resolve_attempted: A flag indicating whether an
                attempt to resolve the RefSeq ID has been made. Defaults to False.
            bgc_path: The path to the downloaded BGC file for
                the genome. Defaults to "".
        """
        self.original_id = original_id
        self.resolved_refseq_id = "" if resolved_refseq_id == "None" else resolved_refseq_id
        self.resolve_attempted = resolve_attempted
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
            "resolved_refseq_id": self.resolved_refseq_id,
            "resolve_attempted": self.resolve_attempted,
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
        # get the best available ID from the dict
        genome_id_data = genome_record["genome_ID"]
        raw_genome_id = get_best_available_genome_id(genome_id_data)
        if raw_genome_id is None or len(raw_genome_id) == 0:
            logger.warning(f'Invalid input genome record "{genome_record}"')
            continue

        # check if genome ID exist in the genome status file
        if raw_genome_id not in gs_dict:
            gs_dict[raw_genome_id] = GenomeStatus(raw_genome_id)

        gs_obj = gs_dict[raw_genome_id]

        logger.info(
            f"Checking for antismash data {i + 1}/{len(genome_records)}, "
            f"current genome ID={raw_genome_id}"
        )
        # first, check if BGC data is downloaded
        if gs_obj.bgc_path and Path(gs_obj.bgc_path).exists():
            logger.info(f"Genome ID {raw_genome_id} already downloaded to {gs_obj.bgc_path}")
            continue
        # second, check if lookup attempted previously
        if gs_obj.resolve_attempted:
            logger.info(f"Genome ID {raw_genome_id} skipped due to previous failed attempt")
            continue

        # if not downloaded or lookup attempted, then try to resolve the ID
        # and download
        logger.info(f"Start lookup process for genome ID {raw_genome_id}")
        gs_obj.resolved_refseq_id = _resolve_refseq_id(genome_id_data)
        gs_obj.resolve_attempted = True

        if gs_obj.resolved_refseq_id == "":
            # give up on this one
            logger.warning(f"Failed lookup for genome ID {raw_genome_id}")
            continue

        # if resolved id is valid, try to download and extract antismash data
        try:
            download_and_extract_antismash_data(
                gs_obj.resolved_refseq_id, project_download_root, project_extract_root
            )

            gs_obj.bgc_path = str(
                Path(project_download_root, gs_obj.resolved_refseq_id + ".zip").absolute()
            )

            output_path = Path(project_extract_root, "antismash", gs_obj.resolved_refseq_id)
            if output_path.exists():
                Path.touch(output_path / "completed", exist_ok=True)

        except Exception:
            gs_obj.bgc_path = ""

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


def _ncbi_genbank_search(genbank_id: str, retry_times: int = 3) -> Tag | NavigableString | None:
    url = NCBI_LOOKUP_URL.format(genbank_id)
    retry = 1
    while retry <= retry_times:
        logger.info(f"Looking up GenBank data for {genbank_id} at {url}")
        resp = httpx.get(url, follow_redirects=True)
        if resp.status_code == httpx.codes.OK:
            # the page should contain a <dl> element with class "assembly_summary_new". retrieving
            # the page seems to fail occasionally in the middle of lengthy sequences of genome
            # lookups, so there might be some throttling going on. this will automatically retry
            # the lookup if the expected content isn't found the first time
            soup = BeautifulSoup(resp.content, "html.parser")
            # find the <dl> element with class "assembly_summary_new"
            dl_element = soup.find("dl", {"class": "assembly_summary_new"})
            if dl_element is not None:
                return dl_element
        retry = retry + 1
        time.sleep(5)

    logger.warning(f"Failed to resolve NCBI genome ID {genbank_id} at URL {url} (after retrying)")
    return None


def _resolve_genbank_accession(genbank_id: str) -> str:
    """Try to get RefSeq id through given GenBank id.

    Args:
        genbank_id: ID for GenBank accession.

    Raises:
        Exception: "Unknown HTML format" if the search of genbank does not give any results.
        Exception: "Expected HTML elements not found" if no match with RefSeq assembly accession is found.

    Returns:
        RefSeq ID if the search is successful, otherwise None.
    """
    logger.info(f"Attempting to resolve Genbank accession {genbank_id} to RefSeq accession")
    # genbank id => genbank seq => refseq

    # The GenBank accession can have several formats:
    # 1: BAFR00000000.1
    # 2: NZ_BAGG00000000.1
    # 3: NC_016887.1
    # Case 1 is the default.
    if "_" in genbank_id:
        # case 2
        if len(genbank_id.split("_")[-1].split(".")[0]) == 12:
            genbank_id = genbank_id.split("_")[-1]
        # case 3
        else:
            genbank_id = genbank_id.lower()

    # get rid of any extraneous whitespace
    genbank_id = genbank_id.strip()
    logger.info(f'Parsed GenBank ID to "{genbank_id}"')

    # run a search using the GenBank accession ID
    try:
        dl_element = _ncbi_genbank_search(genbank_id)
        if dl_element is None or isinstance(dl_element, NavigableString):
            raise Exception("Unknown HTML format")

        refseq_idx = -1
        for field_idx, field in enumerate(dl_element.children):
            # this is the element immediately preceding the one with
            # the actual RefSeq ID we want
            if field.getText().strip() == "RefSeq assembly accession:":
                refseq_idx = field_idx + 1

            # this should be True when we've reached the right element
            if field_idx == refseq_idx:
                refseq_id = field.getText()
                # if it has any spaces, take everything up to first one (some have annotations afterwards)
                if refseq_id.find(" ") != -1:
                    refseq_id = refseq_id[: refseq_id.find(" ")]

                return str(refseq_id)

        if refseq_idx == -1:
            raise Exception("Expected HTML elements not found")
    except Exception as e:
        logger.warning(f"Failed resolving GenBank accession {genbank_id}, error {e}")

    return ""


def _resolve_jgi_accession(jgi_id: str) -> str:
    """Try to get RefSeq id through given JGI id.

    Args:
        jgi_id: JGI_Genome_ID for GenBank accession.

    Returns:
        RefSeq ID if search is successful, otherwise None.
    """
    url = JGI_GENOME_LOOKUP_URL.format(jgi_id)
    logger.info(f"Attempting to resolve JGI_Genome_ID {jgi_id} to GenBank accession via {url}")
    # no User-Agent header produces a 403 Forbidden error on this site...
    try:
        resp = httpx.get(
            url, headers={"User-Agent": USER_AGENT}, timeout=10.0, follow_redirects=True
        )
    except httpx.ReadTimeout:
        logger.warning("Timed out waiting for result of JGI_Genome_ID lookup")
        return ""

    soup = BeautifulSoup(resp.content, "html.parser")
    # find the table entry giving the NCBI assembly accession ID
    link = soup.find("a", href=re.compile("https://www.ncbi.nlm.nih.gov/nuccore/.*"))
    if link is None:
        return ""

    return _resolve_genbank_accession(link.text)


def _resolve_refseq_id(genome_id_data: Mapping[str, str]) -> str:
    """Get the RefSeq ID to which the genome accession is linked.

    Check https://pairedomicsdata.bioinformatics.nl/schema.json.

    Args:
        genome_id_data: dictionary containing information
        for each genome record present.

    Returns:
        RefSeq ID if present, otherwise an empty string.
    """
    if "RefSeq_accession" in genome_id_data:
        # best case, can use this directly
        return genome_id_data["RefSeq_accession"]
    if "GenBank_accession" in genome_id_data:
        # resolve via NCBI
        return _resolve_genbank_accession(genome_id_data["GenBank_accession"])
    if "JGI_Genome_ID" in genome_id_data:
        # resolve via JGI => NCBI
        return _resolve_jgi_accession(genome_id_data["JGI_Genome_ID"])

    logger.warning(f"Unable to resolve genome_ID: {genome_id_data}")
    return ""
