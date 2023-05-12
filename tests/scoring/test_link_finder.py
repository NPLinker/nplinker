import numpy as np
import pandas as pd
from pandas.util.testing import assert_frame_equal
import pytest
from pytest import fixture
from nplinker.scoring.linking import DataLinks
from nplinker.scoring.linking import LinkFinder


@fixture(scope='module')
def linkfinder() -> LinkFinder:
    return LinkFinder()


@fixture(scope='module')
def datalinks(gcfs, spectra, mfs, strains):
    """DataLinks object. See `test_data_links.py` for its actual values."""
    return DataLinks(gcfs, spectra, mfs, strains)


def test_init(linkfinder):
    assert_frame_equal(linkfinder.raw_score_spec_gcf, pd.DataFrame())
    assert_frame_equal(linkfinder.raw_score_mf_gcf, pd.DataFrame())
    assert linkfinder.metcalf_mean is None
    assert linkfinder.metcalf_std is None


def test_cal_score_raw_score(linkfinder, datalinks):
    """Test `cal_score` method for `raw_score_spec_gcf` and `raw_score_mf_gcf`.

    The expected values are calculated manually by using values from `test_init`
    of `test_data_links.py` and the default scoring weights.
    """
    # link type = 'spec-gcf'
    linkfinder.cal_score(datalinks, link_type='spec-gcf')
    assert_frame_equal(
        linkfinder.raw_score_spec_gcf,
        pd.DataFrame([[12, -9, 11], [-9, 12, 11], [1, 1, 21]],
                     index=['spectrum1', 'spectrum2', 'spectrum3'],
                     columns=['gcf1', 'gcf2', 'gcf3']))
    # link type = 'mf-gcf'
    linkfinder.cal_score(datalinks, link_type='mf-gcf')
    assert_frame_equal(
        linkfinder.raw_score_mf_gcf,
        pd.DataFrame([[12, -9, 11], [-9, 12, 11], [1, 1, 21]],
                     index=['mf1', 'mf2', 'mf3'],
                     columns=['gcf1', 'gcf2', 'gcf3']))


def test_cal_score_mean_std(linkfinder, datalinks):
    """Test `cal_score` method for `metcalf_mean` and `metcalf_std`."""
    linkfinder.cal_score(datalinks, link_type='spec-gcf')
    assert isinstance(linkfinder.metcalf_mean, np.ndarray)
    assert isinstance(linkfinder.metcalf_std, np.ndarray)
    assert linkfinder.metcalf_mean.shape == (4, 4
                                             )  # (n_strains+1 , n_strains+1)
    assert linkfinder.metcalf_mean.shape == (4, 4)
    # TODO CG: add tests for values after refactoring _cal_mean_std method
    # assert linkfinder.metcalf_mean == expected_array


def test_get_links_gcf(linkfinder, datalinks, gcfs):
    """Test `get_links` method for input GCF objects."""
    linkfinder.cal_score(datalinks, link_type='spec-gcf')
    linkfinder.cal_score(datalinks, link_type='mf-gcf')
    index_names = ['source', 'target', 'score']

    # cutoff = negative infinity (float)
    links = linkfinder.get_links(*gcfs, score_cutoff=np.NINF)
    assert len(links) == 2
    # expected values got from `test_cal_score_raw_score`
    assert_frame_equal(
        links[0],
        pd.DataFrame([['gcf1', 'gcf2', 'gcf3'] * 3,
                      [
                          *['spectrum1'] * 3,
                          *['spectrum2'] * 3,
                          *['spectrum3'] * 3,
                      ], [12, -9, 11, -9, 12, 11, 1, 1, 21]],
                     index=index_names))
    assert_frame_equal(
        links[1],
        pd.DataFrame([['gcf1', 'gcf2', 'gcf3'] * 3,
                      [
                          *['mf1'] * 3,
                          *['mf2'] * 3,
                          *['mf3'] * 3,
                      ], [12, -9, 11, -9, 12, 11, 1, 1, 21]],
                     index=index_names))

    # cutoff = 0
    links = linkfinder.get_links(*gcfs, score_cutoff=0)
    assert len(links) == 2
    assert_frame_equal(
        links[0],
        pd.DataFrame([['gcf1', 'gcf3', 'gcf2', 'gcf3', 'gcf1', 'gcf2', 'gcf3'],
                      [
                          *['spectrum1'] * 2,
                          *['spectrum2'] * 2,
                          *['spectrum3'] * 3,
                      ], [12, 11, 12, 11, 1, 1, 21]],
                     index=index_names))
    assert_frame_equal(
        links[1],
        pd.DataFrame([['gcf1', 'gcf3', 'gcf2', 'gcf3', 'gcf1', 'gcf2', 'gcf3'],
                      [
                          *['mf1'] * 2,
                          *['mf2'] * 2,
                          *['mf3'] * 3,
                      ], [12, 11, 12, 11, 1, 1, 21]],
                     index=index_names))


def test_get_links_spec(linkfinder, datalinks, spectra):
    """Test `get_links` method for input Spectrum objects."""
    linkfinder.cal_score(datalinks, link_type='spec-gcf')
    linkfinder.cal_score(datalinks, link_type='mf-gcf')
    index_names = ['source', 'target', 'score']
    # cutoff = negative infinity (float)
    links = linkfinder.get_links(*spectra, score_cutoff=np.NINF)
    assert len(links) == 1
    assert_frame_equal(
        links[0],
        pd.DataFrame([[
            *['spectrum1'] * 3,
            *['spectrum2'] * 3,
            *['spectrum3'] * 3,
        ], ['gcf1', 'gcf2', 'gcf3'] * 3, [12, -9, 11, -9, 12, 11, 1, 1, 21]],
                     index=index_names))
    # cutoff = 0
    links = linkfinder.get_links(*spectra, score_cutoff=0)
    assert_frame_equal(
        links[0],
        pd.DataFrame([[
            *['spectrum1'] * 2,
            *['spectrum2'] * 2,
            *['spectrum3'] * 3,
        ], ['gcf1', 'gcf3', 'gcf2', 'gcf3', 'gcf1', 'gcf2', 'gcf3'],
                      [12, 11, 12, 11, 1, 1, 21]],
                     index=index_names))


def test_get_links_mf(linkfinder, datalinks, mfs):
    """Test `get_links` method for input MolecularFamily objects."""
    linkfinder.cal_score(datalinks, link_type='spec-gcf')
    linkfinder.cal_score(datalinks, link_type='mf-gcf')
    index_names = ['source', 'target', 'score']
    # cutoff = negative infinity (float)
    links = linkfinder.get_links(*mfs, score_cutoff=np.NINF)
    assert len(links) == 1
    assert_frame_equal(
        links[0],
        pd.DataFrame([[
            *['mf1'] * 3,
            *['mf2'] * 3,
            *['mf3'] * 3,
        ], ['gcf1', 'gcf2', 'gcf3'] * 3, [12, -9, 11, -9, 12, 11, 1, 1, 21]],
                     index=index_names))
    # cutoff = 0
    links = linkfinder.get_links(*mfs, score_cutoff=0)
    assert_frame_equal(
        links[0],
        pd.DataFrame([[
            *['mf1'] * 2,
            *['mf2'] * 2,
            *['mf3'] * 3,
        ], ['gcf1', 'gcf3', 'gcf2', 'gcf3', 'gcf1', 'gcf2', 'gcf3'],
                      [12, 11, 12, 11, 1, 1, 21]],
                     index=index_names))


@pytest.mark.parametrize("objects, expected", [([], "Empty input objects"),
                                               ("", "Empty input objects")])
def test_get_links_invalid_value(linkfinder, objects, expected):
    with pytest.raises(ValueError) as e:
        linkfinder.get_links(*objects)
    assert expected in str(e.value)


@pytest.mark.parametrize("objects, expected",
                         [([1], "Invalid type {<class 'int'>}"),
                          ([1, 2], "Invalid type {<class 'int'>}"),
                          ("12", "Invalid type {<class 'str'>}")])
def test_get_links_invalid_type(linkfinder, objects, expected):
    with pytest.raises(TypeError) as e:
        linkfinder.get_links(*objects)
    assert expected in str(e.value)


def test_get_links_invalid_mixed_types(linkfinder, spectra, mfs):
    objects = (*spectra, *mfs)
    with pytest.raises(TypeError) as e:
        linkfinder.get_links(*objects)
    assert "Invalid type" in str(e.value)
    assert ".MolecularFamily" in str(e.value)
    assert ".Spectrum" in str(e.value)
