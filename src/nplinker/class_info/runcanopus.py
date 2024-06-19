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
import logging
import os
import subprocess
import sys


logger = logging.getLogger(__name__)


def run_canopus(mgf_file, output_path, extra_params="--maxmz 600 formula zodiac structure canopus"):
    """Runs canopus from the sirius workflow.

    Args:
        mgf_file: str, path to the mgf file with spectra info
        output_path: str, path to canopus_dir for results
        extra_params: str, the extra parameters for running canopus. this
            should always end with canopus, default:
            --maxmz 600 formula zodiac structure canopus
    Returns:
        True if everything runs okay

    Within NPLinker a normal (default) will look like:
    sirius -i METABOLOMICS-SNETS-V2-<id>-download_clustered_spectra-main.mgf
        -o canopus --maxmz 600 formula zodiac structure canopus

    Note that SIRIUS might take a long time when setting a high --maxmz cutoff
    """
    logger.info(
        'run_canopus: input="{}", output="{}", extra_params={}"'.format(
            mgf_file, output_path, extra_params
        )
    )
    if "canopus" not in extra_params and " C " not in extra_params:
        logger.warning("canopus not in sirius command, please check canopus parameters")

    if os.path.exists(os.path.join(output_path, "completed")):
        logger.info("CANOPUS appears to have been run already, skipping!")
        logger.info("To force re-run, delete {}".format(os.path.join(output_path, "completed")))
        return True

    try:
        subprocess.run(["sirius", "-h"], capture_output=True)
    except Exception as e:
        raise Exception(f"Failed to find/run SIRIUS (CANOPUS) (err={e})")

    if not os.path.exists(mgf_file):
        raise Exception(f'mgf input file "{mgf_file}" does not exist!')

    # configure the IO-related parameters
    args = ["sirius", "-i", mgf_file, "-o", output_path]
    if "--log" not in args:
        args.extend(["--log", "INFO"])

    # append the user supplied params, if any
    if len(extra_params) > 0:
        args.extend(extra_params.split(" "))

    logger.info(f"CANOPUS command: {args}")
    result = subprocess.run(args, stdout=sys.stdout, stderr=sys.stderr)
    logger.info("CANOPUS completed with return code {}".format(result.returncode))
    # use subprocess.CompletedProcess.check_returncode() to test if the CANOPUS
    # process exited successfully. This throws an exception for non-zero returncodes
    # which will indicate to the PODPDownloader module that something went wrong.
    result.check_returncode()

    # use presence of this file as a quick way to check if a previous run
    # finished or not
    open(os.path.join(output_path, "completed"), "w").close()

    return True
