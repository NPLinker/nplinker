import sys
import subprocess
import os

from ..logconfig import LogConfig
logger = LogConfig.getLogger(__file__)

# NOTE: for simplicity this is currently written with assumption it will only be
# called in context of nplinker Docker image, where bigscape should be available

def run_bigscape(bigscape_py_path, antismash_path, output_path, pfam_path, cutoffs=[0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7], mibig_14=False, mibig_13=False, extra_params=None):
    logger.info('run_bigscape: input="{}", output="{}", cutoffs={}"'.format(antismash_path, output_path, cutoffs))

    if os.path.exists(os.path.join(output_path, 'completed')):
        logger.info('BiG-SCAPE appears to have been run already, skipping!')
        logger.info('To force re-run, delete {}'.format(os.path.join(output_path, 'completed')))
        return True

    try:
        subprocess.run([bigscape_py_path, '-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except Exception as e:
        raise Exception('Failed to find/run bigscape.py (path={}, err={})'.format(bigscape_py_path, e))

    if not os.path.exists(antismash_path):
        raise Exception('antismash_path "{}" does not exist!'.format(antismash_path))

    if cutoffs is None or len(cutoffs) == 0:
        raise Exception('Should provide at least one cutoff value!')

    if mibig_14 and mibig_13:
        raise Exception('Only one of mibig_14 and mibig_13 should be set!')

    # configure the IO-related parameters, including pfam_dir
    args = [bigscape_py_path, '-i', antismash_path, '-o', output_path, '--pfam_dir', pfam_path]

    # use default cutoffs if not overridden by user
    if extra_params is not None and extra_params.find('cutoffs') == -1:
        args.append('--cutoffs')
        args.append(','.join(map(str, cutoffs)))

    # append the user supplied params, if any
    if extra_params is not None and len(extra_params) > 0:
        args.append(extra_params)

    # TODO mibig 2.0???
    if mibig_13:
        args.append('--mibig13')
    elif mibig_14:
        args.append('--mibig')

    # disable extra clustering (TODO OK?)
    args.append('--clans-off')

    logger.info('bigscape command: {}'.format(args))
    result = subprocess.run(args, stdout=sys.stdout, stderr=sys.stderr)
    logger.info('bigscape completed!')

    # use presence of this file as a quick way to check if a previous run
    # finished or not 
    open(os.path.join(output_path, 'completed'), 'w').close()

    return True

if __name__ == "__main__":
    run_bigscape(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], cutoffs=[0.3])
