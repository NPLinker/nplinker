import os
from pathlib import Path
import pytest
from nplinker.nplinker import NPLinker
from . import DATA_DIR


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


@pytest.mark.skipif(os.environ.get('CI') == 'true',
                    reason="Skip when running on CI")
def test_get_links(npl: NPLinker):
    mc = npl.scoring_method('metcalf')
    mc.cutoff = 3.5
    mc.standardised = True

    actual = npl.get_links(npl.gcfs, mc, and_mode=True)
    assert len(actual) == len(actual.sources) == len(actual.links) == 101

    actual.filter_links(lambda link: link[mc] > 5.0)
    assert len(actual.links) == 60
