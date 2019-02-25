import os
import glob

def _find_via_glob(path, file_type):
    try:
        filename = glob.glob(path)[0]
        return filename
    except OSError as ose:
        raise Exception('{}: unable to find {} in path "{}"'.format(ose, file_type, path))

def setup_paths(root_dir):
    """Given a root directory for a dataset, attempts to construct the filenames
        that are necessary to load everything into nplinker.
    """

    # just <root_dir>/mibig_json
    MIBIG_JSON_DIR = os.path.join(root_dir, 'mibig_json')

    # should have an MGF file named <root_dir>/spectra/<something>.mgf
    spec_dir = os.path.join(root_dir, 'spectra')
    MGF_FILE = _find_via_glob(os.path.join(spec_dir, '*.mgf'), 'MGF file')

    # annotation file(s) may not exist
    ANNOTATION_FILES = glob.glob(os.path.join(spec_dir, '*.annotations.*'))

    # both of these should exist
    EDGES_FILE = _find_via_glob(os.path.join(spec_dir, "*.pairsinfo"), 'edges file')
    NODES_FILE = _find_via_glob(os.path.join(spec_dir, '*.out'), 'nodes file')

    ROOT_PATH = os.path.join(root_dir, 'bigscape')
    ANTISMASH_DIR = os.path.join(root_dir, 'antismash')

    return MIBIG_JSON_DIR, NODES_FILE, ANNOTATION_FILES, MGF_FILE, EDGES_FILE, ROOT_PATH, ANTISMASH_DIR
