# Copyright 2021 The NPLinker Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
import subprocess
import sys
from os import PathLike
from ..logconfig import LogConfig

logger = LogConfig.getLogger(__name__)

# NOTE: for simplicity this is currently written with assumption it will only be
# called in context of nplinker Docker image, where bigscape should be available


def run_bigscape(
    bigscape_py_path: str | PathLike,
    antismash_path: str | PathLike,
    output_path: str | PathLike,
    pfam_path: str | PathLike,
    extra_params: str,
):
    logger.info(
        f'run_bigscape: input="{antismash_path}", output="{output_path}", extra_params={extra_params}"'
    )

    if os.path.exists(os.path.join(output_path, "completed")):
        logger.info("BiG-SCAPE appears to have been run already, skipping!")
        logger.info("To force re-run, delete {%s}", os.path.join(output_path, "completed"))
        return True

    try:
        subprocess.run([bigscape_py_path, "-h"], capture_output=True, check=True)
    except Exception as e:
        raise Exception(f"Failed to find/run bigscape.py (path={bigscape_py_path}, err={e})") from e

    if not os.path.exists(antismash_path):
        raise Exception(f'antismash_path "{antismash_path}" does not exist!')

    # configure the IO-related parameters, including pfam_dir
    args = [bigscape_py_path, "-i", antismash_path, "-o", output_path, "--pfam_dir", pfam_path]

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

    # use presence of this file as a quick way to check if a previous run
    # finished or not
    with open(os.path.join(output_path, "completed"), "w") as f:
        f.close()

    return True


def podp_run_bigscape(
    project_file_cache: str | PathLike,
    PFAM_PATH: str | PathLike,
    do_bigscape: bool,
    extra_bigscape_parameters,
):
    # TODO this currently assumes docker environment, allow customisation?
    # can check if in container with: https://stackoverflow.com/questions/20010199/how-to-determine-if-a-process-runs-inside-lxc-docker
    if not do_bigscape:
        logger.info("BiG-SCAPE disabled by configuration, not running it")
        return

    logger.info('Running BiG-SCAPE! extra_bigscape_parameters="%s"', extra_bigscape_parameters)
    try:
        run_bigscape(
            "bigscape.py",
            os.path.join(project_file_cache, "antismash"),
            os.path.join(project_file_cache, "bigscape"),
            PFAM_PATH,
            extra_bigscape_parameters,
        )
    except Exception as e:
        logger.warning('Failed to run BiG-SCAPE on antismash data, error was "%s"', e)
