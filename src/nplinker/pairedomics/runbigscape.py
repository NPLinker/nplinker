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
from ..logconfig import LogConfig


logger = LogConfig.getLogger(__file__)

# NOTE: for simplicity this is currently written with assumption it will only be
# called in context of nplinker Docker image, where bigscape should be available


def run_bigscape(bigscape_py_path, antismash_path, output_path, pfam_path,
                 extra_params):
    logger.info(
        'run_bigscape: input="{}", output="{}", extra_params={}"'.format(
            antismash_path, output_path, extra_params))

    if os.path.exists(os.path.join(output_path, 'completed')):
        logger.info('BiG-SCAPE appears to have been run already, skipping!')
        logger.info('To force re-run, delete {}'.format(
            os.path.join(output_path, 'completed')))
        return True

    try:
        subprocess.run([bigscape_py_path, '-h'],
                       capture_output=True)
    except Exception as e:
        raise Exception(
            'Failed to find/run bigscape.py (path={}, err={})'.format(
                bigscape_py_path, e))

    if not os.path.exists(antismash_path):
        raise Exception(
            f'antismash_path "{antismash_path}" does not exist!')

    # configure the IO-related parameters, including pfam_dir
    args = [
        bigscape_py_path, '-i', antismash_path, '-o', output_path,
        '--pfam_dir', pfam_path
    ]

    # append the user supplied params, if any
    if len(extra_params) > 0:
        args.extend(extra_params.split(' '))

    logger.info(f'BiG-SCAPE command: {args}')
    result = subprocess.run(args, stdout=sys.stdout, stderr=sys.stderr)
    logger.info('BiG-SCAPE completed with return code {}'.format(
        result.returncode))
    # use subprocess.CompletedProcess.check_returncode() to test if the BiG-SCAPE
    # process exited successfully. This throws an exception for non-zero returncodes
    # which will indicate to the Downloader module that something went wrong.
    result.check_returncode()

    # use presence of this file as a quick way to check if a previous run
    # finished or not
    open(os.path.join(output_path, 'completed'), 'w').close()

    return True


if __name__ == "__main__":
    run_bigscape(sys.argv[1],
                 sys.argv[2],
                 sys.argv[3],
                 sys.argv[4],
                 cutoffs=[0.3])
