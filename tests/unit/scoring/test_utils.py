import pandas as pd
from pandas.testing import assert_frame_equal
from nplinker.scoring.utils import get_presence_gcf_strain
from nplinker.scoring.utils import get_presence_mf_strain
from nplinker.scoring.utils import get_presence_spec_strain
from nplinker.scoring.utils import isinstance_all


def test_isinstance_all():
    assert isinstance_all(1, 2, 3, objtype=int)
    assert not isinstance_all(1, 2, 3, objtype=str)
    assert not isinstance_all(1, 2, "3", objtype=int)


#
# Test get_presence_* functions
#


def test_get_presence_gcf_strain(gcfs, strains):
    presence_gcf_strain = get_presence_gcf_strain(gcfs, strains)
    assert_frame_equal(
        presence_gcf_strain,
        pd.DataFrame(
            [[1, 0, 0], [0, 1, 0], [1, 1, 0]],
            index=["gcf1", "gcf2", "gcf3"],
            columns=["strain1", "strain2", "strain3"],
        ),
    )


def test_get_presence_spec_strain(spectra, strains):
    presence_spec_strain = get_presence_spec_strain(spectra, strains)
    assert_frame_equal(
        presence_spec_strain,
        pd.DataFrame(
            [[1, 0, 0], [0, 1, 0], [1, 1, 0]],
            index=["spectrum1", "spectrum2", "spectrum3"],
            columns=["strain1", "strain2", "strain3"],
        ),
    )


def test_get_presence_mf_strain(mfs, strains):
    presence_mf_strain = get_presence_mf_strain(mfs, strains)
    assert_frame_equal(
        presence_mf_strain,
        pd.DataFrame(
            [[1, 0, 0], [0, 1, 0], [1, 1, 0]],
            index=["mf1", "mf2", "mf3"],
            columns=["strain1", "strain2", "strain3"],
        ),
    )
