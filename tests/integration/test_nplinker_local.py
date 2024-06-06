import hashlib
from pathlib import Path
import pytest
from nplinker.nplinker import NPLinker
from . import DATA_DIR


# Only tests related to data arranging and loading should be put here.
# For tests on scoring/links, add them to `scoring/test_nplinker_scoring.py`.


def get_file_hash(file_path):
    h = hashlib.sha256()
    with open(file_path, "rb") as file:
        while True:
            # Reading is buffered, so we can read smaller chunks.
            chunk = file.read(h.block_size)
            if not chunk:
                break
            h.update(chunk)

    return h.hexdigest()


@pytest.fixture(scope="module")
def npl() -> NPLinker:
    npl = NPLinker(DATA_DIR / "nplinker_local_mode.toml")
    npl.load_data()
    # remove cached score results before running tests
    root_dir = Path(npl.root_dir)
    score_cache = root_dir / "output" / "cache_metcalf_scoring.pckl"
    score_cache.unlink(missing_ok=True)
    return npl


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
    assert len(npl.molfams) == 29
    assert len(npl.strains) == 46
