import os
import json
from os import PathLike
from pathlib import Path
import re
import time
import zipfile
from bs4 import BeautifulSoup
from bs4 import NavigableString
from bs4 import Tag
from deprecated import deprecated
import httpx
from progress.bar import Bar
from nplinker.genomics.antismash import download_and_extract_antismash_data
from nplinker.logconfig import LogConfig

logger = LogConfig.getLogger(__name__)

NCBI_LOOKUP_URL = 'https://www.ncbi.nlm.nih.gov/assembly/?term={}'
JGI_GENOME_LOOKUP_URL = 'https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid={}'
USER_AGENT = 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:86.0) Gecko/20100101 Firefox/86.0'
GENOME_STATUS_FILENAME = "genome_status.json"


class GenomeStatus:
    """A class to represent the status of a single genome.

    The status of genomes is tracked in a JSON file which has a name defined
    in variable `GENOME_STATUS_FILENAME`.
    """

    def __init__(self,
                 original_id: str,
                 resolved_refseq_id: str = "",
                 resolve_attempted: bool = False,
                 bgc_path: str = ""):
        """Initialize a GenomeStatus object for the given genome.

        Args:
            original_id (str): The original ID of the genome.
            resolved_refseq_id (str, optional): The resolved RefSeq ID of the
                genome. Defaults to "".
            resolve_attempted (bool, optional): A flag indicating whether an
                attempt to resolve the RefSeq ID has been made. Defaults to False.
            bgc_path (str, optional): The path to the downloaded BGC file for
                the genome. Defaults to "".
        """
        self.original_id = original_id
        self.resolved_refseq_id = resolved_refseq_id
        self.resolve_attempted = resolve_attempted
        self.bgc_path = bgc_path

    @staticmethod
    def load_from_json(file: str | PathLike) -> dict[str, 'GenomeStatus']:
        """Get a dict of GenomeStatus objects by loading given genome status file.

        Note that an empty dict is returned if the given file doesn't exist.

        Args:
            file(str | PathLike): Path to genome status file.

        Returns:
            dict: dict keys are genome original id and values are GenomeStatus
                objects. An empty dict is returned if the given file doesn't exist.
        """
        genome_status_dict = {}
        if Path(file).exists():
            with open(file, "r") as f:
                data = json.load(f)
            genome_status_dict = {
                gs["original_id"]: GenomeStatus(**gs)
                for gs in data["genome_status"]
            }
        return genome_status_dict

    @staticmethod
    def save_to_json(genome_status_dict: dict[str, 'GenomeStatus'],
                     output_dir: str | PathLike) -> None:
        """Save the given genome status dictionary to a JSON file.

        The JSON file will be saved to the given output directory with the name
        defined in variable `GENOME_STATUS_FILENAME`. The file will be overwritten
        if it already exists.

        Args:
            genome_status_dict (dict[str, 'GenomeStatus']): A dictionary of genome
                status objects to be saved to a JSON file.
            output_dir(str | PathLike): The path to the directory where the JSON
                file will be saved.
        """
        json_data = {
            "genome_status":
            [gs._to_dict() for gs in genome_status_dict.values()],
            "version": "1.0"
        }
        with open(Path(output_dir) / GENOME_STATUS_FILENAME, "w") as f:
            json.dump(json_data, f)

    def _to_dict(self) -> dict:
        """Convert the GenomeStatus object to a dict."""
        return {
            "original_id": self.original_id,
            "resolved_refseq_id": self.resolved_refseq_id,
            "resolve_attempted": self.resolve_attempted,
            "bgc_path": self.bgc_path
        }


def podp_download_and_extract_antismash_data(
        genome_records: list[dict[str, dict[str, str]]],
        project_download_root: str | PathLike,
        project_extract_root: str | PathLike):
    """Download and extract antiSMASH BGC archive for the given genome records.

    Args:
        genome_records(list[dict[str, dict[str, str] | str]]): list of dicts
            representing genome records. The dict of each genome record contains
                - key(str): "genome_ID"
                - value(dict[str, str]): a dict containing information about genome
                type, label and accession ids (RefSeq, GenBank, and/or JGI).
        project_download_root(str | PathLike): Path to the directory to place
            downloaded archive in.
        project_extract_root(str | PathLike): Path to the directory downloaded archive
            will be extracted to.
            Note that an `antismash` directory will be created in the specified
            `extract_root` if it doesn't exist. The files will be extracted to
            `<extract_root>/antismash/<antismash_id>` directory.
    """

    if not Path(project_download_root).exists():
        # otherwise in case of failed first download, the folder doesn't exist and
        # genome_status_file can't be written
        Path(project_download_root).mkdir(parents=True, exist_ok=True)

    genome_status_file = Path(project_download_root, GENOME_STATUS_FILENAME)
    genome_status = _get_genome_status_log(genome_status_file)

    for i, genome_record in enumerate(genome_records):
        # get the best available ID from the dict
        raw_genome_id = _get_best_available_genome_id(
            genome_record['genome_ID'])
        if raw_genome_id is None or len(raw_genome_id) == 0:
            logger.warning(
                f'Ignoring genome record "{genome_record}" due to missing genome ID field'
            )
            continue

        # use this to check if the lookup has already been attempted and if
        # so if the file is cached locally
        if raw_genome_id not in genome_status:
            genome_status[raw_genome_id] = GenomeStatus(raw_genome_id, "None")

        genome_obj = genome_status[raw_genome_id]

        logger.info(
            f'Checking for antismash data {i + 1}/{len(genome_records)}, current genome ID={raw_genome_id}'
        )
        # first check if file is cached locally
        if (genome_obj.filename and Path(genome_obj.filename).exists()):
            # file already downloaded
            logger.info(
                f'Genome ID {raw_genome_id} already downloaded to {genome_obj.filename}'
            )
            continue
        if genome_obj.attempted:
            # lookup attempted previously but failed
            logger.info(
                f'Genome ID {raw_genome_id} skipped due to previous failure')
            continue
        # if no existing file and no lookup attempted, can start process of
        # trying to retrieve the data

        # lookup the ID
        logger.info(f'Beginning lookup process for genome ID {raw_genome_id}')

        genome_obj.resolved_refseq_id = _resolve_refseq_id(
            genome_record['genome_ID'])

        if not isinstance(genome_obj.resolved_refseq_id, str):
            raise TypeError(
                f"genome_obj.resolved_refseq_id should be a string. Instead got: {type(genome_obj.resolved_refseq_id)}"
            )

        genome_obj.attempted = True

        if genome_obj.resolved_refseq_id == "":
            # give up on this one
            logger.warning(f'Failed lookup for genome ID {raw_genome_id}')
            continue

        # if we got a refseq ID, now try to download and extract the data from antismash
        download_and_extract_antismash_data(genome_obj.resolved_refseq_id,
                                            project_download_root,
                                            project_extract_root)

        genome_obj.filename = str(
            Path(project_download_root,
                 genome_obj.resolved_refseq_id + '.zip').absolute())
        output_path = Path(project_extract_root, 'antismash',
                           genome_obj.resolved_refseq_id)
        Path.touch(output_path / 'completed', exist_ok=True)

    missing = len([x for x in genome_status.values() if len(x.filename) == 0])
    logger.info(
        f'Dataset has {missing} missing sets of antiSMASH data (from a total of {len(genome_records)})'
    )

    for obj in genome_status.values():
        obj.to_csv(genome_status_file)

    if missing == len(genome_records):
        logger.warning('Failed to successfully retrieve ANY genome data!')


def _get_genome_status_log(
        genome_status_file: PathLike) -> dict[str, GenomeStatus]:
    """Get a dict of GenomeStatus objects by reading given genome status file.
    Note that a empty dict is returned if the given file does not exist.

    Args:
        genome_status_file(PathLike): Path to genome status file that records
            genome IDs and local filenames to avoid repeating time-consuming
            HTTP requests each time the app is loaded.

    Returns:
        dict: dict keys are genome original id and values are GenomeStatus objects.
    """

    genome_status = {}

    # GENOME_STATUS_FILENAME is read, then in the for loop over the genome records it gets updated,
    # and finally it is saved again in GENOME_STATUS_FILENAME which is overwritten
    if Path(genome_status_file).exists():
        with open(genome_status_file) as f:
            for line in csv.reader(f):
                asobj = GenomeStatus(*line)
                genome_status[asobj.original_id] = asobj

    return genome_status


def _get_best_available_genome_id(
        genome_id_data: dict[str, str]) -> str | None:
    """Get the best available ID from genome_id_data dict.

    Args:
        genome_id_data(dict): dictionary containing information
        for each genome record present.

    Returns:
        str | None: ID for the genome, if present, otherwise None.
    """

    if 'RefSeq_accession' in genome_id_data:
        return genome_id_data['RefSeq_accession']
    if 'GenBank_accession' in genome_id_data:
        return genome_id_data['GenBank_accession']
    if 'JGI_Genome_ID' in genome_id_data:
        return genome_id_data['JGI_Genome_ID']

    logger.warning(
        f'No known genome ID field in genome data: {genome_id_data}')
    return None


def _ncbi_genbank_search(genbank_id: str,
                         retry_times: int = 3) -> Tag | NavigableString | None:

    url = NCBI_LOOKUP_URL.format(genbank_id)
    retry = 1
    while retry <= retry_times:
        logger.debug(f'Looking up GenBank data for {genbank_id} at {url}')
        resp = httpx.get(url, follow_redirects=True)
        if resp.status_code == httpx.codes.OK:
            # the page should contain a <dl> element with class "assembly_summary_new". retrieving
            # the page seems to fail occasionally in the middle of lengthy sequences of genome
            # lookups, so there might be some throttling going on. this will automatically retry
            # the lookup if the expected content isn't found the first time
            soup = BeautifulSoup(resp.content, 'html.parser')
            # find the <dl> element with class "assembly_summary_new"
            dl_element = soup.find('dl', {'class': 'assembly_summary_new'})
            if dl_element is not None:
                return dl_element
        retry = retry + 1
        time.sleep(5)

    logger.warning(
        f'Failed to resolve NCBI genome ID {genbank_id} at URL {url} (after retrying)'
    )
    return None


def _resolve_genbank_accession(genbank_id: str) -> str:
    """Try to get RefSeq id through given GenBank id.

    Args:
        genbank_id(str): ID for GenBank accession.

    Raises:
        Exception: "Unknown HTML format" if the search of genbank does not give any results.
        Exception: "Expected HTML elements not found" if no match with RefSeq assembly accession is found.

    Returns:
        str | None: RefSeq ID if the search is successful, otherwise None.
    """
    logger.info(
        f'Attempting to resolve Genbank accession {genbank_id} to RefSeq accession'
    )
    # genbank id => genbank seq => refseq

    # The GenBank accession can have several formats:
    # 1: BAFR00000000.1
    # 2: NZ_BAGG00000000.1
    # 3: NC_016887.1
    # Case 1 is the default.
    if '_' in genbank_id:
        # case 2
        if len(genbank_id.split('_')[-1].split('.')[0]) == 12:
            genbank_id = genbank_id.split('_')[-1]
        # case 3
        else:
            genbank_id = genbank_id.lower()

    # get rid of any extraneous whitespace
    genbank_id = genbank_id.strip()
    logger.debug(f'Parsed GenBank ID to "{genbank_id}"')

    # run a search using the GenBank accession ID
    try:
        dl_element = _ncbi_genbank_search(genbank_id)
        if dl_element is None:
            raise Exception('Unknown HTML format')

        refseq_idx = -1
        for field_idx, field in enumerate(dl_element.children):
            # this is the element immediately preceding the one with
            # the actual RefSeq ID we want
            if field.getText().strip() == 'RefSeq assembly accession:':
                refseq_idx = field_idx + 1

            # this should be True when we've reached the right element
            if field_idx == refseq_idx:
                refseq_id = field.getText()
                # if it has any spaces, take everything up to first one (some have annotations afterwards)
                if refseq_id.find(' ') != -1:
                    refseq_id = refseq_id[:refseq_id.find(' ')]

                return str(refseq_id)

        if refseq_idx == -1:
            raise Exception('Expected HTML elements not found')
    except Exception as e:
        logger.warning(
            f'Failed resolving GenBank accession {genbank_id}, error {e}')

    return ""


def _resolve_jgi_accession(jgi_id: str) -> str:
    """Try to get RefSeq id through given JGI id.

    Args:
        jgi_id(str): JGI_Genome_ID for GenBank accession.

    Returns:
        str | None: Return RefSeq ID if search is successful, otherwise None.
    """
    url = JGI_GENOME_LOOKUP_URL.format(jgi_id)
    logger.info(
        f'Attempting to resolve JGI_Genome_ID {jgi_id} to GenBank accession via {url}'
    )
    # no User-Agent header produces a 403 Forbidden error on this site...
    try:
        resp = httpx.get(url,
                         headers={'User-Agent': USER_AGENT},
                         timeout=10.0,
                         follow_redirects=True)
    except httpx.ReadTimeout:
        logger.warning('Timed out waiting for result of JGI_Genome_ID lookup')
        return ""

    soup = BeautifulSoup(resp.content, 'html.parser')
    # find the table entry giving the NCBI assembly accession ID
    link = soup.find(
        'a', href=re.compile('https://www.ncbi.nlm.nih.gov/nuccore/.*'))
    if link is None:
        return ""

    return _resolve_genbank_accession(link.text)


def _resolve_refseq_id(genome_id_data: dict[str, str]) -> str:
    """Get the RefSeq ID to which the genome accession is linked.
    Check https://pairedomicsdata.bioinformatics.nl/schema.json.

    Args:
        genome_id_data(dict): dictionary containing information
        for each genome record present.

    Returns:
        str: Return RefSeq ID if present, otherwise an empty string.
    """
    if 'RefSeq_accession' in genome_id_data:
        # best case, can use this directly
        return genome_id_data['RefSeq_accession']
    if 'GenBank_accession' in genome_id_data:
        # resolve via NCBI
        return _resolve_genbank_accession(genome_id_data['GenBank_accession'])
    if 'JGI_Genome_ID' in genome_id_data:
        # resolve via JGI => NCBI
        return _resolve_jgi_accession(genome_id_data['JGI_Genome_ID'])

    logger.warning(f'Unable to resolve genome_ID: {genome_id_data}')
    return ""
