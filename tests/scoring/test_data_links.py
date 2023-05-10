import pandas as pd
from pandas.util.testing import assert_frame_equal
from pytest import fixture
from nplinker.scoring.linking import DataLinks


@fixture
def datalinks(gcfs, spectra, mfs, strains):
    return DataLinks(gcfs, spectra, mfs, strains)


def test_init(datalinks):
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
