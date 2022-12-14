import pytest
from nplinker.nplinker import NPLinker
import os

from . import DATA_DIR

@pytest.fixture(scope='module')
def instance() -> NPLinker:
    npl = NPLinker(str(DATA_DIR / 'nplinker_demo1.toml'))
    npl.load_data()
    return npl

def runningInDocker():
    with open('/proc/self/cgroup', 'r') as procfile:
        for line in procfile:
            fields = line.strip().split('/')
            if 'docker' in fields:
                return True
    return False

@pytest.mark.skipif(os.environ.get('CI') == 'true', reason="Skip when running on CI")
def test_load_data(instance: NPLinker):
    assert len(instance.bgcs) == 390
    assert len(instance.gcfs) == 113
    assert len(instance.spectra) == 25935
    assert len(instance.molfams) == 25769


@pytest.mark.skipif(os.environ.get('CI') == 'true' , reason="Skip when running on CI")
def test_get_links(instance: NPLinker):
    mc = instance.scoring_method('metcalf')
    mc.cutoff = 3.5
    mc.standardised = True

    actual = instance.get_links(instance.gcfs, mc, and_mode=True)
    assert len(actual) == len(actual.sources) == len(actual.links) == 101
    
    actual.filter_links(lambda link: link[mc] > 5.0)
    assert len(actual.links) == 60

