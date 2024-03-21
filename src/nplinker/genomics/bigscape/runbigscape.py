from __future__ import annotations
import os
import subprocess
import sys
from os import PathLike
from ...logconfig import LogConfig


logger = LogConfig.getLogger(__name__)

PFAM_PATH = os.path.join(sys.prefix, "nplinker_lib")


def run_bigscape(
    antismash_path: str | PathLike,
    output_path: str | PathLike,
    extra_params: str,
):
    bigscape_py_path = "bigscape.py"
    logger.info(
        f'run_bigscape: input="{antismash_path}", output="{output_path}", extra_params={extra_params}"'
    )

    try:
        subprocess.run([bigscape_py_path, "-h"], capture_output=True, check=True)
    except Exception as e:
        raise Exception(f"Failed to find/run bigscape.py (path={bigscape_py_path}, err={e})") from e

    if not os.path.exists(antismash_path):
        raise Exception(f'antismash_path "{antismash_path}" does not exist!')

    # configure the IO-related parameters, including pfam_dir
    args = [bigscape_py_path, "-i", antismash_path, "-o", output_path, "--pfam_dir", PFAM_PATH]

    # append the user supplied params, if any
    if len(extra_params) > 0:
        args.extend(extra_params.split(" "))

    logger.info(f"BiG-SCAPE command: {args}")
    result = subprocess.run(args, stdout=sys.stdout, stderr=sys.stderr, check=True)
    logger.info(f"BiG-SCAPE completed with return code {result.returncode}")
    # use subprocess.CompletedProcess.check_returncode() to test if the BiG-SCAPE
    # process exited successfully. This throws an exception for non-zero returncodes
    # which will indicate to the PODPDownloader module that something went wrong.
    result.check_returncode()

    return True
