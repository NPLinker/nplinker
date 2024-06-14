import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal
from nplinker.scoring import MetcalfScoring


def test_init(npl):
    mc = MetcalfScoring()
    assert mc.name == "metcalf"
    assert mc.npl is None
    assert mc.metcalf_weights == (10, -10, 0, 1)
    assert_frame_equal(mc.presence_gcf_strain, pd.DataFrame())
    assert_frame_equal(mc.presence_spec_strain, pd.DataFrame())
    assert_frame_equal(mc.presence_mf_strain, pd.DataFrame())
    assert_frame_equal(mc.raw_score_spec_gcf, pd.DataFrame(columns=["spec", "gcf", "score"]))
    assert_frame_equal(mc.raw_score_mf_gcf, pd.DataFrame(columns=["mf", "gcf", "score"]))
    assert mc.metcalf_mean is None
    assert mc.metcalf_std is None


#
# Test the `setup` method
#


def test_setup(mc, gcfs, spectra, mfs, strains):
    """Test `setup` method when cache file does not exist."""
    assert_frame_equal(
        mc.presence_gcf_strain,
        pd.DataFrame([[1, 0, 0], [0, 1, 0], [1, 1, 0]], index=gcfs, columns=list(strains)),
    )
    assert_frame_equal(
        mc.presence_spec_strain,
        pd.DataFrame(
            [[1, 0, 0], [0, 1, 0], [1, 1, 0]],
            index=spectra,
            columns=list(strains),
        ),
    )
    assert_frame_equal(
        mc.presence_mf_strain,
        pd.DataFrame(
            [[1, 0, 0], [0, 1, 0], [1, 1, 0]],
            index=mfs,
            columns=list(strains),
        ),
    )

    df = pd.DataFrame([[12, -9, 11], [-9, 12, 11], [1, 1, 21]], index=spectra, columns=gcfs)
    df_melted = df.reset_index().melt(id_vars="index")
    df_melted.columns = ["spec", "gcf", "score"]
    assert_frame_equal(mc.raw_score_spec_gcf, df_melted)

    df = pd.DataFrame([[12, -9, 11], [-9, 12, 11], [1, 1, 21]], index=mfs, columns=gcfs)
    df_melted = df.reset_index().melt(id_vars="index")
    df_melted.columns = ["mf", "gcf", "score"]
    assert_frame_equal(mc.raw_score_mf_gcf, df_melted)

    assert isinstance(mc.metcalf_mean, np.ndarray)
    assert isinstance(mc.metcalf_std, np.ndarray)
    assert mc.metcalf_mean.shape == (4, 4)  # (n_strains+1 , n_strains+1)
    assert mc.metcalf_std.shape == (4, 4)


#
# Test the `get_links` method
#


def test_get_links_default(mc, gcfs, spectra, mfs):
    # same as cutoff=0, standardised=False
    lg = mc.get_links()
    assert lg[gcfs[0]][spectra[0]][mc.name].value == 12
    assert lg[gcfs[1]].get(spectra[0]) is None
    assert lg[gcfs[2]][spectra[0]][mc.name].value == 11
    assert lg[gcfs[0]][mfs[0]][mc.name].value == 12
    assert lg[gcfs[1]][mfs[1]][mc.name].value == 12
    assert lg[gcfs[2]][mfs[2]][mc.name].value == 21


@pytest.mark.parametrize(
    "objects, expected",
    [
        ([1], "Invalid type <class 'int'>. .*"),
        ([1, 2], "Invalid type <class 'int'>. .*"),
        ("12", "Invalid type <class 'str'>. .*"),
    ],
)
def test_get_links_invalid_input_type(mc, objects, expected):
    with pytest.raises(TypeError, match=expected):
        mc.get_links(*objects)


def test_get_links_invalid_mixed_types(mc, spectra, mfs):
    objects = (*spectra, *mfs)
    with pytest.raises(TypeError, match="Input objects must be of the same type."):
        mc.get_links(*objects)


def test_get_links_gcf_standardised_false(mc, gcfs, spectra, mfs):
    """Test `get_links` method when input is GCF objects and `standardised` is False."""
    # when cutoff is negative infinity, i.e. taking all scores
    lg = mc.get_links(*gcfs, cutoff=np.NINF, standardised=False)
    assert lg[gcfs[0]][spectra[0]][mc.name].value == 12
    assert lg[gcfs[1]][spectra[0]][mc.name].value == -9
    assert lg[gcfs[2]][spectra[0]][mc.name].value == 11
    assert lg[gcfs[0]][mfs[0]][mc.name].value == 12
    assert lg[gcfs[1]][mfs[1]][mc.name].value == 12
    assert lg[gcfs[2]][mfs[2]][mc.name].value == 21

    # when test cutoff is 0, i.e. taking scores >= 0
    lg = mc.get_links(*gcfs, cutoff=0, standardised=False)
    assert lg[gcfs[0]][spectra[0]][mc.name].value == 12
    assert lg[gcfs[1]].get(spectra[0]) is None
    assert lg[gcfs[2]][spectra[0]][mc.name].value == 11
    assert lg[gcfs[0]][mfs[0]][mc.name].value == 12
    assert lg[gcfs[1]][mfs[1]][mc.name].value == 12
    assert lg[gcfs[2]][mfs[2]][mc.name].value == 21


def test_get_links_gcf_standardised_true(mc, gcfs):
    """Test `get_links` method when input is GCF objects and `standardised` is True."""
    lg = mc.get_links(*gcfs, cutoff=np.NINF, standardised=True)
    assert len(lg.links) == 18

    lg = mc.get_links(*gcfs, cutoff=0, standardised=True)
    assert len(lg.links) == 14


def test_get_links_spec_standardised_false(mc, gcfs, spectra):
    """Test `get_links` method when input is Spectrum objects and `standardised` is False."""
    lg = mc.get_links(*spectra, cutoff=np.NINF, standardised=False)
    assert lg[spectra[0]][gcfs[0]][mc.name].value == 12
    assert lg[spectra[0]][gcfs[1]][mc.name].value == -9
    assert lg[spectra[0]][gcfs[2]][mc.name].value == 11

    lg = mc.get_links(*spectra, cutoff=0, standardised=False)
    assert lg[spectra[0]][gcfs[0]][mc.name].value == 12
    assert lg[spectra[0]].get(gcfs[1]) is None
    assert lg[spectra[0]][gcfs[2]][mc.name].value == 11


def test_get_links_spec_standardised_true(mc, gcfs, spectra):
    """Test `get_links` method when input is Spectrum objects and `standardised` is True."""
    lg = mc.get_links(*spectra, cutoff=np.NINF, standardised=True)
    assert len(lg.links) == 9

    lg = mc.get_links(*spectra, cutoff=0, standardised=True)
    assert len(lg.links) == 7


def test_get_links_mf_standardised_false(mc, gcfs, mfs):
    """Test `get_links` method when input is MolecularFamily objects and `standardised` is False."""
    lg = mc.get_links(*mfs, cutoff=np.NINF, standardised=False)
    assert lg[mfs[0]][gcfs[0]][mc.name].value == 12
    assert lg[mfs[0]][gcfs[1]][mc.name].value == -9
    assert lg[mfs[0]][gcfs[2]][mc.name].value == 11

    lg = mc.get_links(*mfs, cutoff=0, standardised=False)
    assert lg[mfs[0]][gcfs[0]][mc.name].value == 12
    assert lg[mfs[0]].get(gcfs[1]) is None
    assert lg[mfs[0]][gcfs[2]][mc.name].value == 11


def test_get_links_mf_standardised_true(mc, gcfs, mfs):
    """Test `get_links` method when input is MolecularFamily objects and `standardised` is True."""
    lg = mc.get_links(*mfs, cutoff=np.NINF, standardised=True)
    assert len(lg.links) == 9

    lg = mc.get_links(*mfs, cutoff=0, standardised=True)
    assert len(lg.links) == 7
