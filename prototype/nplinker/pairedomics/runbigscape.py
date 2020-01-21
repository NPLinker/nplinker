import sys
import subprocess
import os

# NOTE: for simplicity this is currently written with assumption it will only be
# called in context of nplinker Docker image, where bigscape should be available

def run_bigscape(antismash_path, output_path, cutoffs=[0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7], mibig_14=False, mibig_13=False):
    try:
        subprocess.run(['bigscape.py', '-h'])
    except Exception as e:
        raise Exception('Failed to find bigscape.py ({})'.format(e))

    if not os.path.exists(antismash_path):
        raise Exception('antismash_path "{}" does not exist!'.format(antismash_path))

    if len(cutoffs) == 0:
        raise Exception('Should provide at least one cutoff value!')

    if mibig_14 and mibig_13:
        raise Exception('Only one of mibig_14 and mibig_13 should be set!')

    args = ['bigscape.py', '-i', antismash_path, '-o', output_path, '--cutoffs', ','.join(map(str, cutoffs))]
    if mibig_13:
        args.append('--mibig13')
    if mibig_14:
        args.append('--mibig')
    # TODO mibig 2.0???


    print(args) 
    result = subprocess.run(args, stdout=sys.stdout, stderr=sys.stderr)

    print('BiG-SCAPE completed running!')
    open(os.path.join(output_path, 'completed'), 'w').close()

    return True

if __name__ == "__main__":
    run_bigscape()
