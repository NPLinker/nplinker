import os
from pathlib import Path
import pytest
from nplinker.nplinker import NPLinker
from . import DATA_DIR


# NOTE: This file only contains tests that run locally and are skipped on CI.
# Basically, only tests related to data loading should be put here.
# For tests on scoring/links, add them to `scoring/test_nplinker_scoring.py`.


@pytest.fixture(scope='module')
def npl() -> NPLinker:
    npl = NPLinker(str(DATA_DIR / 'nplinker_demo1.toml'))
    npl.load_data()
    # remove cached results before running tests
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
