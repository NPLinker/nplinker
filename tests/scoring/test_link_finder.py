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
    return DataLinks(gcfs, spectra, mfs, strains)


def test_init(linkfinder):
    assert_frame_equal(linkfinder.raw_score_spec_gcf, pd.DataFrame())
    assert_frame_equal(linkfinder.raw_score_mf_gcf, pd.DataFrame())
    assert linkfinder.metcalf_mean is None
    assert linkfinder.metcalf_std is None


def test_cal_score(linkfinder, datalinks):
    linkfinder.cal_score(datalinks, link_type='spec-gcf')
    assert_frame_equal(
        linkfinder.raw_score_spec_gcf,
        pd.DataFrame([[12, -9, 11], [-9, 12, 11], [1, 1, 21]],
                     index=['spectrum1', 'spectrum2', 'spectrum3'],
                     columns=['gcf1', 'gcf2', 'gcf3']))

    linkfinder.cal_score(datalinks, link_type='mf-gcf')
    assert_frame_equal(
        linkfinder.raw_score_mf_gcf,
        pd.DataFrame([[12, -9, 11], [-9, 12, 11], [1, 1, 21]],
                     index=['mf1', 'mf2', 'mf3'],
                     columns=['gcf1', 'gcf2', 'gcf3']))

    # TODO CG: add tests for values after refactoring _cal_mean_std
    assert isinstance(linkfinder.metcalf_mean, np.ndarray)
    assert isinstance(linkfinder.metcalf_std, np.ndarray)
    assert linkfinder.metcalf_mean.shape == (4, 4
                                             )  # (n_strains+1 , n_strains+1)
    assert linkfinder.metcalf_mean.shape == (4, 4)
    # generated metcalf_mean
    # array([[  3.,   2.,   1.,   0.],
    #    [ -8.,  -2.,   4.,  10.],
    #    [-19.,  -6.,   7.,  20.],
    #    [-30., -10.,  10.,  30.]])
    # generated metcalf_std
    # array([[1.        , 1.        , 1.        , 1.        ],
    #    [1.        , 9.89949494, 9.89949494, 1.        ],
    #    [1.        , 9.89949494, 9.89949494, 1.        ],
    #    [1.        , 1.        , 1.        , 1.        ]])


def test_get_links_gcf(linkfinder, datalinks, gcfs):
    linkfinder.cal_score(datalinks, link_type='spec-gcf')
    linkfinder.cal_score(datalinks, link_type='mf-gcf')
    index_names = ['source', 'target', 'score']
    # cutoff = negative infinity (float)
    links = linkfinder.get_links(*gcfs, score_cutoff=np.NINF)
    assert len(links) == 2
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


def test_get_links_exceptions(linkfinder):
    with pytest.raises(TypeError) as e:
        linkfinder.get_links("")
    assert "Invalid type {<class 'str'>}" in str(e.value)
