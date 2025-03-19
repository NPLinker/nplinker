from __future__ import annotations
import logging
from os import PathLike
from pathlib import Path
from typing import Optional
import requests


logger = logging.getLogger(__name__)


def submit_antismash_job(genbank_filepath: str | PathLike) -> str:
    """Submits an antiSMASH job using the provided GenBank file.

    This function sends a GenBank file to the antiSMASH API
    and retrieves the job ID if the submission is successful.

    Args:
        genbank_filepath (str | PathLike): The path to the GenBank file to be submitted.

    Returns:
        str: The job ID if the submission.

    Raises:
        requests.exceptions.RequestException: If there is an issue with the HTTP request.
        RuntimeError: If the API response does not contain a job ID.
    """
    url = "https://antismash.secondarymetabolites.org/api/v1.0/submit"
    genbank_filepath = Path(genbank_filepath)

    with open(genbank_filepath, "rb") as file:
        files = {"seq": file}
        response = requests.post(url, files=files)
        response.raise_for_status()  # Raise an exception for HTTP errors

    data = response.json()
    if "id" not in data:
        raise RuntimeError("No antiSMASH job ID returned")
    return data["id"]


def query_antismash_job(job_id: str) -> Optional[dict]:
    """Gets the status of an antiSMASH job.

    Args:
        job_id (str): The job ID to query.

    Returns:
        dict: The response JSON if successful, otherwise None.
    """
    url = f"https://antismash.secondarymetabolites.org/api/v1.0/status/{job_id}"

    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()  # Raise an exception for HTTP errors
        return response.json()

    except requests.exceptions.RequestException as req_err:
        logger.error(f"Request failed for job_id {job_id}: {req_err}")
    except ValueError as json_err:  # Handles JSON decoding errors
        logger.error(f"Invalid JSON response for job_id {job_id}: {json_err}")
    except Exception as err:
        logger.error(f"Unexpected error while getting job state for job_id {job_id}: {err}")


def antismash_job_is_done(job_id: str) -> bool:
    """Checks if the antiSMASH job is complete by polling the job status.

    Args:
        job_id (str): The job ID to query.

    Returns:
        bool: True if the job is done, False if the job is still running.

    Raises:
        RuntimeError: If the job status could not be retrieved or if the job failed.
        ValueError: If the job state is missing or unexpected in the response.
    """
    response = query_antismash_job(job_id)

    if response is None:
        raise RuntimeError(f"Failed to retrieve job status for job_id {job_id}")
    if "state" not in response:
        raise ValueError(f"Job state missing in response for job_id: {job_id}")

    job_state = response["state"]
    if job_state in ("running", "queued"):
        return False
    if job_state == "done":
        return True
    if job_state == "failed":
        job_status = response.get("status", "No error message provided")
        raise RuntimeError(f"AntiSMASH job {job_id} failed with an error: {job_status}")
    else:
        raise ValueError(
            f"Unexpected job state for antismash job ID {job_id}. Job state: {job_state}"
        )
