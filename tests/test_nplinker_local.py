import hashlib
import os
from pathlib import Path
import pytest
from nplinker.nplinker import NPLinker
from . import DATA_DIR


# NOTE: This file only contains tests that run locally and are skipped on CI.
# Basically, only tests related to data loading should be put here.
# For tests on scoring/links, add them to `scoring/test_nplinker_scoring.py`.


def get_file_hash(file_path):
    h = hashlib.sha256()
    with open(file_path, 'rb') as file:
        while True:
            # Reading is buffered, so we can read smaller chunks.
            chunk = file.read(h.block_size)
            if not chunk:
                break
            h.update(chunk)

    return h.hexdigest()


@pytest.fixture(scope='module')
def npl() -> NPLinker:
    npl = NPLinker(str(DATA_DIR / 'nplinker_demo1.toml'))
    npl.load_data()
    hash_proj_file = get_file_hash(
        os.path.join(Path(npl._loader._root).parent.parent,
                     npl._loader._platform_id + '.json'))
    if hash_proj_file != '22e4f20d6f8aa425b2040479d0b6c00e7d3deb03f8fc4a277b3b91eb07c9ad72':
        pytest.exit(
            'PoDP project file has changed, please clean your local cache folder and rerun the tests.'
        )
    # remove cached score results before running tests
    root_dir = Path(npl.root_dir)
    score_cache = root_dir / 'metcalf' / 'metcalf_scores.pckl'
    score_cache.unlink(missing_ok=True)
    return npl


@pytest.mark.skipif(os.environ.get('CI') == 'true',
                    reason="Skip when running on CI")
def test_load_data(npl: NPLinker):

    assert len(npl.bgcs) == 390
    assert len(npl.gcfs) == 113
    assert len(npl.spectra) == 25935
    assert len(npl.molfams) == 25769
