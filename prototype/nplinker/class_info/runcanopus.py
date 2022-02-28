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


import sys
import subprocess
import os

from ..logconfig import LogConfig
logger = LogConfig.getLogger(__file__)


def run_canopus(sirius_path, mgf_file, output_path,
                extra_params='--maxmz 850 formula zodiac structure canopus'):
    """Runs canopus from the sirius workflow

    Args:
        sirius_path: str, where to find sirius
        mgf_file: str, path to the mgf file with spectra info
        output_path: str, path to canopus_dir for results
        extra_params: str, the extra parameters for running canopus. this
            should always end with canopus, default:
            --maxmz 850 formula zodiac structure canopus
    Returns:
        True if everything runs okay

    Within NPLinker a normal (default) will look like:
    sirius -i METABOLOMICS-SNETS-V2-<id>-download_clustered_spectra-main.mgf
        -o canopus --maxmz 850 formula zodiac structure canopus
    """
    logger.info('run_canopus: input="{}", output="{}", extra_params={}"'.format(mgf_file, output_path, extra_params))
    if 'canopus' not in extra_params or ' C ' not in extra_params:
        logger.warn('canopus not in sirius command, please check canopus parameters')

    if os.path.exists(os.path.join(output_path, 'canopus_summary.tsv')):
        logger.info('CANOPUS appears to have been run already, skipping!')
        logger.info('To force re-run, delete {}'.format(os.path.join(output_path, 'canopus_summary.tsv')))
        return True

    try:
        subprocess.run([sirius_path, '-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except Exception as e:
        raise Exception('Failed to find/run SIRIUS (CANOPUS) (path={}, err={})'.format(sirius_path, e))

    if not os.path.exists(mgf_file):
        raise Exception('mgf input file "{}" does not exist!'.format(mgf_file))

    # configure the IO-related parameters
    args = [sirius_path, '-i', mgf_file, '-o', output_path]

    # append the user supplied params, if any
    if len(extra_params) > 0:
        args.extend(extra_params.split(' '))

    logger.info('CANOPUS command: {}'.format(args))
    result = subprocess.run(args, stdout=sys.stdout, stderr=sys.stderr)
    logger.info('CANOPUS completed with return code {}'.format(result.returncode))
    # use subprocess.CompletedProcess.check_returncode() to test if the CANOPUS
    # process exited successfully. This throws an exception for non-zero returncodes
    # which will indicate to the Downloader module that something went wrong.
    result.check_returncode()

    # use presence of this file as a quick way to check if a previous run
    # finished or not
    open(os.path.join(output_path, 'canopus_summary.tsv'), 'w').close()

    return True


if __name__ == "__main__":
    run_canopus(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
