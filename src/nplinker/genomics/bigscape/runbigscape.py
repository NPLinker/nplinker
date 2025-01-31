from __future__ import annotations
import logging
import os
import subprocess
import sys
from os import PathLike
from typing import Literal


logger = logging.getLogger(__name__)

PFAM_PATH = os.path.join(sys.prefix, "nplinker_lib")


def run_bigscape(
    antismash_path: str | PathLike,
    output_path: str | PathLike,
    extra_params: str,
    version: Literal[1, 2] = 1,
) -> bool:
    """Runs BiG-SCAPE to cluster BGCs.

    The behavior of this function is slightly different depending on the version of
    BiG-SCAPE that is set to run using the configuration file.
    Mostly this means a different set of parameters is used between the two versions.

    The AntiSMASH output directory should be a directory that contains GBK files.
    The directory can contain subdirectories, in which case BiG-SCAPE will search
    recursively for GBK files. E.g.:

    ```
    example_folder
        ├── organism_1
        │  ├── organism_1.region001.gbk
        │  ├── organism_1.region002.gbk
        │  ├── organism_1.region003.gbk
        │  ├── organism_1.final.gbk          <- skipped!
        │  └── ...
        ├── organism_2
        │  ├── ...
        └── ...
    ```

    By default, only GBK Files with "cluster" or "region" in the filename are
    accepted. GBK Files with "final" in the filename are excluded.

    Args:
        antismash_path: Path to the antismash output directory.
        output_path: Path to the output directory where BiG-SCAPE will write its results.
        extra_params: Additional parameters to pass to BiG-SCAPE.
        version: The version of BiG-SCAPE to run. Must be 1 or 2.

    Returns:
        True if BiG-SCAPE ran successfully, False otherwise.

    Raises:
        ValueError: If an unexpected BiG-SCAPE version number is specified.
        FileNotFoundError: If the antismash_path does not exist or if the BiG-SCAPE python
            script could not be found.
        RuntimeError: If BiG-SCAPE fails to run.

    Examples:
        >>>  from nplinker.genomics.bigscape import run_bigscape
        >>> run_bigscape(antismash_path="./antismash", output_path="./output",
        ... extra_params="--help", version=1)
    """
    # switch to correct version of BiG-SCAPE
    if version == 1:
        bigscape_py_path = "bigscape.py"
    elif version == 2:
        bigscape_py_path = "bigscape-v2.py"
    else:
        raise ValueError("Invalid BiG-SCAPE version number. Expected: 1 or 2.")

    try:
        subprocess.run([bigscape_py_path, "-h"], capture_output=True, check=True)
    except Exception as e:
        raise FileNotFoundError(
            f"Failed to find/run BiG-SCAPE executable program (path={bigscape_py_path}, err={e})"
        ) from e

    if not os.path.exists(antismash_path):
        raise FileNotFoundError(f'antismash_path "{antismash_path}" does not exist!')

    logger.info(f"Running BiG-SCAPE version {version}")
    logger.info(
        f'run_bigscape: input="{antismash_path}", output="{output_path}", extra_params={extra_params}"'
    )

    # assemble arguments. first argument is the python file
    args = [bigscape_py_path]

    # version 2 points to specific Pfam file, version 1 points to directory
    # version 2 also requires the cluster subcommand
    if version == 1:
        args.extend(["--pfam_dir", PFAM_PATH])
    elif version == 2:
        args.extend(["cluster", "--pfam_path", os.path.join(PFAM_PATH, "Pfam-A.hmm")])

    # add input and output paths. these are unchanged
    args.extend(["-i", str(antismash_path), "-o", str(output_path)])

    # append the user supplied params, if any
    if len(extra_params) > 0:
        args.extend(extra_params.split(" "))

    logger.info(f"BiG-SCAPE command: {args}")
    result = subprocess.run(args, stdout=sys.stdout, stderr=sys.stderr)

    # return true on any non-error return code
    if result.returncode == 0:
        logger.info(f"BiG-SCAPE completed with return code {result.returncode}")
        return True

    # otherwise log details and raise a runtime error
    logger.error(f"BiG-SCAPE failed with return code {result.returncode}")
    logger.error(f"output: {str(result.stdout)}")
    logger.error(f"stderr: {str(result.stderr)}")

    raise RuntimeError(f"Failed to run BiG-SCAPE with error code {result.returncode}")
