import logging
import re
from typing import Any
from typing import Callable
from typing import Literal
from typing import Mapping
import httpx
from bs4 import BeautifulSoup


JGI_GENOME_LOOKUP_URL = (
    "https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid={}"
)
USER_AGENT = "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:86.0) Gecko/20100101 Firefox/86.0"

logger = logging.getLogger(__name__)


def get_latest_assembly_accession(acc: str) -> str:
    """Retrieve the latest NCBI genome assembly accession for a given accession.

    This function retrieves the most recent genome assembly accession
    from the revision history of the provided accession. It prioritizes
    RefSeq accessions, and if unavailable, falls back to GenBank accessions.

    Args:
        acc (str): The accession identifier to resolve.

    Returns:
        str: The latest valid genome assembly accession.

    Raises:
        ValueError: If no valid RefSeq or GenBank accession is found in
            the assembly revision history.
    """
    revision_history = _get_revision_history(acc)

    acc_priority = ("refseq_accession", "genbank_accession")
    for acc_type in acc_priority:
        assembly_revisions = revision_history.get("assembly_revisions", [])
        revisions_with_acc = [entry for entry in assembly_revisions if acc_type in entry]
        if revisions_with_acc:
            latest_revision = max(revisions_with_acc, key=lambda x: x["release_date"])
            return str(latest_revision[acc_type])

    raise ValueError("No valid genome accession found in assembly revision history")


def resolve_genome_accession(genome_id_data: Mapping[Any, Any]) -> str:
    """Gets the NCBI genome assembly accession.

    Gets the latest RefSeq genome assembly accession, or if not available,
    the latest GenBank genome assembly accession.

    This function gets the latest RefSeq genome assembly accession, or if not
    available, the latest GenBank genome assembly accession. It attempts to get
    a genome accession by checking the genome id date for specific ID types in
    the following order:
    1. RefSeq_accession
    2. GenBank_accession
    3. JGI_Genome_ID

    For each ID type, it uses a corresponding resolver function to process
    the ID. If a resolver fails, a warning is logged, and the function
    proceeds to the next ID type. If no valid genome assembly accession can be
    retrieved, a RuntimeError is raised.

    Args:
        genome_id_data (dict): A dictionary containing genome ID types as keys
            and their corresponding values.

    Returns:
        str: The retrieved genome accession, prioritizing RefSeq if available,
            otherwise GenBank.

    Raises:
        RuntimeError: If no valid assembly accessions can be retrieved.

    Logs:
        Warning messages if a resolver fails for a specific ID type.
    """
    resolver_priority = ["RefSeq_accession", "GenBank_accession", "JGI_Genome_ID"]
    resolvers: dict[str, Callable] = {
        "RefSeq_accession": _resolve_refseq,
        "GenBank_accession": _resolve_genbank,
        "JGI_Genome_ID": _resolve_jgi,
    }

    for id_type in resolver_priority:
        if id_type not in genome_id_data:
            continue

        resolver = resolvers[id_type]
        try:
            genome_id = genome_id_data[id_type].strip()
            return str(resolver(genome_id))
        except Exception as e:
            logger.warning(f"Failed to resolve {id_type}: {e}")

    raise RuntimeError("No valid assembly accessions found")


def validate_assembly_acc(acc: str, acc_type: Literal["RefSeq", "GenBank"]) -> None:
    """Validates a NCBI genome assembly accession string based on its type.

    This function checks if the provided genome assembly accession string
    adheres to the expected format for the specified accession type
    ('RefSeq' or 'GenBank'). It ensures that the accession starts with the
    correct prefix ('GCF_' for RefSeq or 'GCA_' for GenBank) and contains
    a version number.

    Args:
        acc (str): The genome assembly accession string to validate.
        acc_type (Literal["RefSeq", "GenBank"]): The type of accession, either
            'RefSeq' or 'GenBank'.

    Raises:
        ValueError: If the accession type is invalid.
        ValueError: If the accession does not start with the expected prefix.
        ValueError: If the accession is missing a version number.
    """
    if acc_type == "RefSeq":
        prefix = "GCF_"
    elif acc_type == "GenBank":
        prefix = "GCA_"
    else:
        raise ValueError(
            f"Invalid genome assembly accession type '{acc_type}'. Expected 'RefSeq' or 'GenBank'."
        )
    if not acc.startswith(prefix):
        raise ValueError(f"Invalid {acc_type} assembly accession (must start with {prefix}): {acc}")
    if "." not in acc:
        raise ValueError(f"Invalid assembly accession (missing version number): {acc}")


def _get_revision_history(assembly_acc: str) -> dict[str, Any]:
    """Fetches the revision history of a genome assembly from the NCBI Datasets API.

    Args:
        assembly_acc (str): The accession number of the genome assembly.

    Returns:
        dict[str, Any]: A dictionary containing the revision history of the specified assembly.

    Raises:
        httpx.HTTPStatusError: If the HTTP request fails or returns a non-success status code.
        ValueError: If no revision history data is found for the given accession number.
    """
    url = (
        f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{assembly_acc}/revision_history"
    )
    resp = httpx.get(url, headers={"User-Agent": USER_AGENT}, timeout=10.0, follow_redirects=True)
    resp.raise_for_status()
    revision_history = resp.json()
    if not revision_history:
        raise ValueError(f"No Assembly Revision data found for {assembly_acc}")
    if not isinstance(revision_history, dict):
        raise ValueError(f"Unexpected response format: {type(revision_history)}")
    return revision_history


def _resolve_refseq(acc: str) -> str:
    """Retrieves the latest NCBI genome assembly accession for a given RefSeq accession.

    This function validates the provided RefSeq accession and retrieves the
    most up-to-date RefSeq genome assembly accession, or if not found,
    alternatively the most up-to-date RefSeq GenBank genome assembly accession.

    Args:
        acc (str): The RefSeq accession to resolve.

    Returns:
        str: The latest NCBI assembly accession corresponding to the given RefSeq accession.

    Raises:
        ValueError: If the provided accession is invalid or not a RefSeq accession.
    """
    validate_assembly_acc(acc, "RefSeq")
    return get_latest_assembly_accession(acc)


def _resolve_genbank(acc: str) -> str:
    """Retrieves the latest NCBI genome assembly accession for a given GenBank accession.

    This function validates the provided GenBank accession and retrieves
    most up-to-date RefSeq genome assembly accession, or if not found,
    alternatively the most up-to-date RefSeq GenBank genome assembly accession.

    Args:
        acc (str): The GenBank accession to resolve.

    Returns:
        str: The latest NCBI assembly accession corresponding to the given GenBank accession.

    Raises:
        ValueError: If the provided accession is invalid or not a RefSeq accession.
    """
    validate_assembly_acc(acc, "GenBank")
    return get_latest_assembly_accession(acc)


def _resolve_jgi(jgi_genome_id: str) -> str:
    """Resolves a JGI Genome ID to its corresponding NCBI assembly accession.

    This function queries a predefined JGI genome lookup URL using the provided
    JGI Genome ID. It parses the HTML response to extract the NCBI assembly
    accession link. If no valid link is found, an exception is raised. The
    function then retrieves the most up-to-date NCBI assembly accession for the
    found assembly accession.

    Args:
        jgi_genome_id (str): The JGI Genome ID to resolve.

    Returns:
        str: The latest NCBI assembly accession corresponding to the given JGI Genome ID.

    Raises:
        ValueError: If no NCBI accessions can be found for the given JGI Genome ID.
        httpx.HTTPStatusError: If the HTTP request to the JGI genome lookup URL fails.
    """
    url = JGI_GENOME_LOOKUP_URL.format(jgi_genome_id)
    resp = httpx.get(url, headers={"User-Agent": USER_AGENT}, timeout=10.0, follow_redirects=True)
    resp.raise_for_status()

    soup = BeautifulSoup(resp.content, "html.parser")
    link = soup.find("a", href=re.compile("https://www.ncbi.nlm.nih.gov/datasets/genome/.*"))
    if not link:
        raise ValueError(f"Unable to find NCBI accessions for the JGI Genome ID: {jgi_genome_id}. ")
    assembly_acc = link.text
    return get_latest_assembly_accession(assembly_acc)
