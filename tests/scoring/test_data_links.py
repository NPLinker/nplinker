import pandas as pd
from pandas.util.testing import assert_frame_equal
import pytest
from pytest import fixture
from nplinker.metabolomics.singleton_family import SingletonFamily
from nplinker.scoring.linking import DataLinks


@fixture(scope='module')
def datalinks(gcfs, spectra, mfs, strains) -> DataLinks:
    return DataLinks(gcfs, spectra, mfs, strains)


def test_init(datalinks):
    """Test that the DataLinks object is initialised correctly.

        Multiple private methods are called in the init method, so we test
        that the correct dataframes are created.
    """
    # test occorrences
    col_names = ['strain1', 'strain2', 'strain3']
    assert_frame_equal(
        datalinks.occurrence_gcf_strain,
        pd.DataFrame([[1, 0, 0], [0, 1, 0], [1, 1, 0]],
                     index=['gcf1', 'gcf2', 'gcf3'],
                     columns=col_names))
    assert_frame_equal(
        datalinks.occurrence_spec_strain,
        pd.DataFrame([[1, 0, 0], [0, 1, 0], [1, 1, 0]],
                     index=['spectrum1', 'spectrum2', 'spectrum3'],
                     columns=col_names))
    assert_frame_equal(
        datalinks.occurrence_mf_strain,
        pd.DataFrame([[1, 0, 0], [0, 1, 0], [1, 1, 0]],
                     index=['mf1', 'mf2', 'mf3'],
                     columns=col_names))
    # test co-occorrences spec-gcf
    col_names = ['gcf1', 'gcf2', 'gcf3']
    assert_frame_equal(
        datalinks.cooccurrence_spec_gcf,
        pd.DataFrame([[1, 0, 1], [0, 1, 1], [1, 1, 2]],
                     index=['spectrum1', 'spectrum2', 'spectrum3'],
                     columns=col_names))
    assert_frame_equal(
        datalinks.cooccurrence_spec_notgcf,
        pd.DataFrame([[0, 1, 0], [1, 0, 0], [1, 1, 0]],
                     index=['spectrum1', 'spectrum2', 'spectrum3'],
                     columns=col_names))
    assert_frame_equal(
        datalinks.cooccurrence_notspec_gcf,
        pd.DataFrame([[0, 1, 1], [1, 0, 1], [0, 0, 0]],
                     index=['spectrum1', 'spectrum2', 'spectrum3'],
                     columns=col_names))
    assert_frame_equal(
        datalinks.cooccurrence_notspec_notgcf,
        pd.DataFrame([[2, 1, 1], [1, 2, 1], [1, 1, 1]],
                     index=['spectrum1', 'spectrum2', 'spectrum3'],
                     columns=col_names))
    # test co-occorrences mf-gcf
    assert_frame_equal(
        datalinks.cooccurrence_mf_gcf,
        pd.DataFrame([[1, 0, 1], [0, 1, 1], [1, 1, 2]],
                     index=['mf1', 'mf2', 'mf3'],
                     columns=col_names))
    assert_frame_equal(
        datalinks.cooccurrence_mf_notgcf,
        pd.DataFrame([[0, 1, 0], [1, 0, 0], [1, 1, 0]],
                     index=['mf1', 'mf2', 'mf3'],
                     columns=col_names))
    assert_frame_equal(
        datalinks.cooccurrence_notmf_gcf,
        pd.DataFrame([[0, 1, 1], [1, 0, 1], [0, 0, 0]],
                     index=['mf1', 'mf2', 'mf3'],
                     columns=col_names))
    assert_frame_equal(
        datalinks.cooccurrence_notmf_notgcf,
        pd.DataFrame([[2, 1, 1], [1, 2, 1], [1, 1, 1]],
                     index=['mf1', 'mf2', 'mf3'],
                     columns=col_names))


def test_get_common_strains_spec(datalinks, spectra, gcfs, strains_list):
    """Test get_common_strains method for input spectra and gcfs."""
    sut = datalinks.get_common_strains(spectra[:2], gcfs)
    expected = {
        (spectra[0], gcfs[0]): [strains_list[0]],
        (spectra[0], gcfs[1]): [],
        (spectra[0], gcfs[2]): [strains_list[0]],
        (spectra[1], gcfs[0]): [],
        (spectra[1], gcfs[1]): [strains_list[1]],
        (spectra[1], gcfs[2]): [strains_list[1]]
    }
    assert sut == expected

    sut = datalinks.get_common_strains(spectra[:2],
                                       gcfs,
                                       filter_no_shared=True)
    expected = {
        (spectra[0], gcfs[0]): [strains_list[0]],
        (spectra[0], gcfs[2]): [strains_list[0]],
        (spectra[1], gcfs[1]): [strains_list[1]],
        (spectra[1], gcfs[2]): [strains_list[1]]
    }
    assert sut == expected


def test_get_common_strains_mf(datalinks, mfs, gcfs, strains_list):
    """Test get_common_strains method for input mfs and gcfs."""
    sut = datalinks.get_common_strains(mfs[:2], gcfs)
    expected = {
        (mfs[0], gcfs[0]): [strains_list[0]],
        (mfs[0], gcfs[1]): [],
        (mfs[0], gcfs[2]): [strains_list[0]],
        (mfs[1], gcfs[0]): [],
        (mfs[1], gcfs[1]): [strains_list[1]],
        (mfs[1], gcfs[2]): [strains_list[1]]
    }
    assert sut == expected

    sut = datalinks.get_common_strains(mfs[:2], gcfs, filter_no_shared=True)
    expected = {
        (mfs[0], gcfs[0]): [strains_list[0]],
        (mfs[0], gcfs[2]): [strains_list[0]],
        (mfs[1], gcfs[1]): [strains_list[1]],
        (mfs[1], gcfs[2]): [strains_list[1]]
    }
    assert sut == expected


def test_get_common_strains_sf(datalinks, mfs, gcfs, strains_list):
    """Test get_common_strains method for input SingletonFamily."""
    smf = SingletonFamily()

    sut = datalinks.get_common_strains([smf], gcfs)
    assert sut == {}

    # the expected are same as `test_get_common_strains_mf`
    mfs_mix = (*mfs[:2], smf)
    sut = datalinks.get_common_strains(mfs_mix, gcfs)
    expected = {
        (mfs[0], gcfs[0]): [strains_list[0]],
        (mfs[0], gcfs[1]): [],
        (mfs[0], gcfs[2]): [strains_list[0]],
        (mfs[1], gcfs[0]): [],
        (mfs[1], gcfs[1]): [strains_list[1]],
        (mfs[1], gcfs[2]): [strains_list[1]]
    }
    assert sut == expected

    sut = datalinks.get_common_strains(mfs_mix, gcfs, filter_no_shared=True)
    expected = {
        (mfs[0], gcfs[0]): [strains_list[0]],
        (mfs[0], gcfs[2]): [strains_list[0]],
        (mfs[1], gcfs[1]): [strains_list[1]],
        (mfs[1], gcfs[2]): [strains_list[1]]
    }
    assert sut == expected


@pytest.mark.parametrize(
    "first_arg, expected",
    [([], "Empty list for first or second argument."),
     ((), "Empty list for first or second argument."),
     ([1], "First argument must be Spectrum or MolecularFamily objects.")])
def test_get_common_strains_exception_first_arg(datalinks, spectra, mfs, gcfs,
                                                first_arg, expected):
    """Test get_common_strains method for invalid 1st arugment."""
    with pytest.raises(ValueError) as e:
        datalinks.get_common_strains(first_arg, gcfs)
    assert expected in str(e.value)

    # mixed input
    with pytest.raises(ValueError) as e:
        datalinks.get_common_strains(spectra + mfs, gcfs)
    assert "First argument must be Spectrum or MolecularFamily objects." in str(
        e.value)


@pytest.mark.parametrize("second_arg, expected",
                         [([], "Empty list for first or second argument."),
                          ((), "Empty list for first or second argument."),
                          ([1], "Second argument must be GCF objects.")])
def test_get_common_strains_exception_second_arg(datalinks, spectra,
                                                 second_arg, expected):
    """Test get_common_strains method for invalid 2nd argument."""
    with pytest.raises(ValueError) as e:
        datalinks.get_common_strains(spectra, second_arg)
    assert expected in str(e.value)
