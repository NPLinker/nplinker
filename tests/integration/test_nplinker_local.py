import os
import pickle
import pytest
from nplinker.genomics import GCF
from nplinker.metabolomics import MolecularFamily
from nplinker.metabolomics import Spectrum
from nplinker.nplinker import NPLinker
from . import DATA_DIR


@pytest.fixture(scope="module")
def npl(root_dir) -> NPLinker:
    os.environ["NPLINKER_ROOT_DIR"] = root_dir
    npl = NPLinker(DATA_DIR / "nplinker_local_mode.toml")
    npl.load_data()
    return npl


def test_init(npl, root_dir):
    assert str(npl.config.root_dir) == root_dir
    assert npl.config.mode == "local"
    assert npl.config.log.level == "DEBUG"

    assert npl.root_dir == root_dir


# ---------------------------------------------------------------------------------------------------
# After manually checking data files for PODP MSV000079284, we have the following numbers:
# 370 BGCs from antismash files
# 114 GCFs, including:
#   - 49 singleton GCFs
#   - 1 mibig-only GCF (not singleton)
#   - 12 GCFs (neither singleton nor mibig-only) have mibig bgcs and in total 20 mibig BGCs are used
# 25935 spectra, including:
#   - 24652 spectra have strain info (from strain mapping file)
#   - 1283 spectra do not have strain info
# 25769 molecular families, including:
#   - 25740 singleton families
#   - 29 non-singleton families
# 26 strains from strain mapping file
# ---------------------------------------------------------------------------------------------------
# So, after data loading, we should get following numbers in the tests:
# 390 BGCs = 370 antismash BGCs + 20 mibig BGCs
# 64 GCFs (neither singleton nor mibig-only) = 114 GCFs - 49 singleton GCFs - 1 mibig-only GCF
# 24652 spectra (having strain info)
# 29 molecular families (non-singleton)
# 46 strains = 26 strains from strain mapping file + 20 strains from mibig
# ---------------------------------------------------------------------------------------------------


def test_load_data(npl: NPLinker):
    assert len(npl.bgcs) == 390
    assert len(npl.gcfs) == 64
    assert len(npl.spectra) == 24652
    assert len(npl.mfs) == 29
    assert len(npl.strains) == 46


def test_get_links(npl):
    # default scoring parameters are used (cutoff=0, standardised=False),
    # so all score values should be >= 0
    scoring_method = "metcalf"
    lg = npl.get_links(npl.gcfs[:3], scoring_method)
    for _, _, scores in lg.links:
        score = scores[scoring_method]
        assert score.value >= 0

    lg = npl.get_links(npl.spectra[:1], scoring_method)
    for _, _, scores in lg.links:
        score = scores[scoring_method]
        assert score.value >= 0

    lg = npl.get_links(npl.mfs[:1], scoring_method)
    for _, _, scores in lg.links:
        score = scores[scoring_method]
        assert score.value >= 0


def test_save_data(npl):
    scoring_method = "metcalf"
    links = npl.get_links(npl.gcfs[:3], scoring_method)

    pickle_file = os.path.join(npl.output_dir, "npl.pkl")
    npl.save_data(pickle_file, links)

    with open(pickle_file, "rb") as f:
        bgcs, gcfs, spectra, mfs, strains, lg = pickle.load(f)

    # tests from `test_load_data`
    assert len(bgcs) == 390
    assert len(gcfs) == 64
    assert len(spectra) == 24652
    assert len(mfs) == 29
    assert len(strains) == 46

    # tests from `test_get_links`
    for obj1, obj2, scores in lg.links:
        score = scores[scoring_method]
        assert score.value >= 0

        if isinstance(obj1, GCF):
            assert obj1 in gcfs
        elif isinstance(obj1, Spectrum):
            assert obj1 in spectra
        elif isinstance(obj1, MolecularFamily):
            assert obj1 in mfs
        else:
            assert False
