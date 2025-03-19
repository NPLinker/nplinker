from __future__ import annotations
import logging
from os import PathLike
from pathlib import Path
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
    return str(data["id"])


def antismash_job_is_done(job_id: str) -> bool:
    """Determines if the antiSMASH job has completed by checking its status.

    This function queries the antiSMASH API to retrieve the current state
    of the job and determines whether it has finished successfully, is still
    in progress, or has encountered an error.

    Args:
        job_id (str): The unique identifier of the antiSMASH job.

    Returns:
        bool: True if the job is completed successfully, False if it is still
            running or queued.

    Raises:
        RuntimeError: If the job has failed or if the API response indicates an error.
        ValueError: If the job state is missing or an unexpected state is encountered
            in the API response.
        requests.exceptions.HTTPError: If an HTTP error occurs during the API request.
    """
    url = f"https://antismash.secondarymetabolites.org/api/v1.0/status/{job_id}"

    response = requests.get(url, timeout=10)
    response.raise_for_status()  # Raise exception for HTTP errors
    respose_data = response.json()

    if "state" not in respose_data:
        raise ValueError(f"Job state missing in response for job_id: {job_id}")

    job_state = respose_data["state"]
    if job_state in ("running", "queued"):
        return False
    if job_state == "done":
        return True
    if job_state == "failed":
        job_status = respose_data.get("status", "No error message provided")
        raise RuntimeError(f"AntiSMASH job {job_id} failed with an error: {job_status}")
    else:
        raise ValueError(
            f"Unexpected job state for antismash job ID {job_id}. Job state: {job_state}"
        )
