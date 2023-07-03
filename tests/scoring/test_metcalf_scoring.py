import numpy as np
from numpy.testing import assert_array_equal
from pandas.util.testing import assert_frame_equal
import pytest
from nplinker.scoring import LinkCollection
from nplinker.scoring import MetcalfScoring
from nplinker.scoring import ObjectLink
from nplinker.scoring.linking import DataLinks
from nplinker.scoring.linking import LinkFinder


def test_init(npl):
    mc = MetcalfScoring(npl)
    assert mc.npl == npl
    assert mc.name == 'metcalf'
    assert mc.cutoff == 1.0
    assert mc.standardised is True
    assert mc.DATALINKS is None
    assert mc.LINKFINDER is None


def test_setup(mc, datalinks, linkfinder):
    """Test `setup` method when cache file does not exist."""
    assert isinstance(mc.DATALINKS, DataLinks)
    assert isinstance(mc.LINKFINDER, LinkFinder)

    assert_frame_equal(mc.DATALINKS.occurrence_gcf_strain,
                       datalinks.occurrence_gcf_strain)
    assert_frame_equal(mc.DATALINKS.cooccurrence_spec_gcf,
                       datalinks.cooccurrence_spec_gcf)

    assert_frame_equal(mc.LINKFINDER.raw_score_spec_gcf,
                       linkfinder.raw_score_spec_gcf)
    assert_frame_equal(mc.LINKFINDER.raw_score_mf_gcf,
                       linkfinder.raw_score_mf_gcf)
    assert_array_equal(mc.LINKFINDER.metcalf_mean, linkfinder.metcalf_mean)
    assert_array_equal(mc.LINKFINDER.metcalf_std, linkfinder.metcalf_std)


def test_setup_load_cache(mc, npl, datalinks, linkfinder, caplog):
    """Test `setup` method when cache file exists."""
    mc.setup(npl)
    assert "MetcalfScoring.setup loading cached data" in caplog.text
    assert "MetcalfScoring.setup caching results" not in caplog.text

    assert isinstance(mc.DATALINKS, DataLinks)
    assert isinstance(mc.LINKFINDER, LinkFinder)

    assert_frame_equal(mc.DATALINKS.occurrence_gcf_strain,
                       datalinks.occurrence_gcf_strain)
    assert_frame_equal(mc.DATALINKS.cooccurrence_spec_gcf,
                       datalinks.cooccurrence_spec_gcf)

    assert_frame_equal(mc.LINKFINDER.raw_score_spec_gcf,
                       linkfinder.raw_score_spec_gcf)
    assert_frame_equal(mc.LINKFINDER.raw_score_mf_gcf,
                       linkfinder.raw_score_mf_gcf)
    assert_array_equal(mc.LINKFINDER.metcalf_mean, linkfinder.metcalf_mean)
    assert_array_equal(mc.LINKFINDER.metcalf_std, linkfinder.metcalf_std)


def test_get_links_gcf_standardised_false(mc, gcfs, spectra, mfs):
    """Test `get_links` method when input is GCF objects and `standardised` is False."""
    # test raw scores (no standardisation)
    mc.standardised = False

    # when cutoff is negative infinity, i.e. taking all scores
    mc.cutoff = np.NINF
    links = mc.get_links(*gcfs, link_collection=LinkCollection())
    assert isinstance(links, LinkCollection)
    links = links.links  # dict of link values
    assert len(links) == 3
    assert {i.gcf_id for i in links.keys()} == {'gcf1', 'gcf2', 'gcf3'}
    assert isinstance(links[gcfs[0]][spectra[0]], ObjectLink)
    # expected values are from `test_get_links_gcf` of test_link_finder.py
    assert links[gcfs[0]][spectra[0]].data(mc) == 12
    assert links[gcfs[1]][spectra[0]].data(mc) == -9
    assert links[gcfs[2]][spectra[0]].data(mc) == 11
    assert links[gcfs[0]][mfs[0]].data(mc) == 12
    assert links[gcfs[1]][mfs[1]].data(mc) == 12
    assert links[gcfs[2]][mfs[2]].data(mc) == 21

    # when test cutoff is 0, i.e. taking scores >= 0
    mc.cutoff = 0
    links = mc.get_links(*gcfs, link_collection=LinkCollection())
    assert isinstance(links, LinkCollection)
    links = links.links
    assert {i.gcf_id for i in links.keys()} == {'gcf1', 'gcf2', 'gcf3'}
    assert isinstance(links[gcfs[0]][spectra[0]], ObjectLink)
    assert links[gcfs[0]][spectra[0]].data(mc) == 12
    assert links[gcfs[1]].get(spectra[0]) is None
    assert links[gcfs[2]][spectra[0]].data(mc) == 11
    assert links[gcfs[0]][mfs[0]].data(mc) == 12
    assert links[gcfs[1]][mfs[1]].data(mc) == 12
    assert links[gcfs[2]][mfs[2]].data(mc) == 21


@pytest.mark.skip(reason='To add after refactoring relevant code.')
def test_get_links_gcf_standardised_true(mc, gcfs, spectra, mfs):
    """Test `get_links` method when input is GCF objects and `standardised` is True."""
    mc.standardised = True
    ...


def test_get_links_spec_standardised_false(mc, gcfs, spectra):
    """Test `get_links` method when input is Spectrum objects and `standardised` is False."""
    mc.standardised = False

    mc.cutoff = np.NINF
    links = mc.get_links(*spectra, link_collection=LinkCollection())
    assert isinstance(links, LinkCollection)
    links = links.links  # dict of link values
    assert len(links) == 3
    assert {i.spectrum_id
            for i in links.keys()} == {'spectrum1', 'spectrum2', 'spectrum3'}
    assert isinstance(links[spectra[0]][gcfs[0]], ObjectLink)
    assert links[spectra[0]][gcfs[0]].data(mc) == 12
    assert links[spectra[0]][gcfs[1]].data(mc) == -9
    assert links[spectra[0]][gcfs[2]].data(mc) == 11

    mc.cutoff = 0
    links = mc.get_links(*spectra, link_collection=LinkCollection())
    assert isinstance(links, LinkCollection)
    links = links.links  # dict of link values
    assert len(links) == 3
    assert {i.spectrum_id
            for i in links.keys()} == {'spectrum1', 'spectrum2', 'spectrum3'}
    assert isinstance(links[spectra[0]][gcfs[0]], ObjectLink)
    assert links[spectra[0]][gcfs[0]].data(mc) == 12
    assert links[spectra[0]].get(gcfs[1]) is None
    assert links[spectra[0]][gcfs[2]].data(mc) == 11


@pytest.mark.skip(reason='To add after refactoring relevant code.')
def test_get_links_spec_standardised_true(mc, gcfs, spectra):
    """Test `get_links` method when input is Spectrum objects and `standardised` is True."""
    mc.standardised = True
    ...


def test_get_links_mf_standardised_false(mc, gcfs, mfs):
    """Test `get_links` method when input is MolecularFamily objects and `standardised` is False."""
    mc.standardised = False

    mc.cutoff = np.NINF
    links = mc.get_links(*mfs, link_collection=LinkCollection())
    assert isinstance(links, LinkCollection)
    links = links.links
    assert len(links) == 3
    assert {i.family_id for i in links.keys()} == {'mf1', 'mf2', 'mf3'}
    assert isinstance(links[mfs[0]][gcfs[0]], ObjectLink)
    assert links[mfs[0]][gcfs[0]].data(mc) == 12
    assert links[mfs[0]][gcfs[1]].data(mc) == -9
    assert links[mfs[0]][gcfs[2]].data(mc) == 11

    mc.cutoff = 0
    links = mc.get_links(*mfs, link_collection=LinkCollection())
    assert isinstance(links, LinkCollection)
    links = links.links
    assert len(links) == 3
    assert {i.family_id for i in links.keys()} == {'mf1', 'mf2', 'mf3'}
    assert isinstance(links[mfs[0]][gcfs[0]], ObjectLink)
    assert links[mfs[0]][gcfs[0]].data(mc) == 12
    assert links[mfs[0]].get(gcfs[1]) is None
    assert links[mfs[0]][gcfs[2]].data(mc) == 11


@pytest.mark.skip(reason='To add after refactoring relevant code.')
def test_get_links_mf_standardised_true(mc, gcfs, mfs):
    """Test `get_links` method when input is MolecularFamily objects and `standardised` is True."""
    mc.standardised = True
    ...


@pytest.mark.parametrize("objects, expected", [([], "Empty input objects"),
                                               ("", "Empty input objects")])
def test_get_links_invalid_input_value(mc, objects, expected):
    with pytest.raises(ValueError) as e:
        mc.get_links(*objects, link_collection=LinkCollection())
    assert expected in str(e.value)


@pytest.mark.parametrize("objects, expected",
                         [([1], "Invalid type {<class 'int'>}"),
                          ([1, 2], "Invalid type {<class 'int'>}"),
                          ("12", "Invalid type {<class 'str'>}")])
def test_get_links_invalid_input_type(mc, objects, expected):
    with pytest.raises(TypeError) as e:
        mc.get_links(*objects, link_collection=LinkCollection())
    assert expected in str(e.value)


def test_get_links_invalid_mixed_types(mc, spectra, mfs):
    objects = (*spectra, *mfs)
    with pytest.raises(TypeError) as e:
        mc.get_links(*objects, link_collection=LinkCollection())
    assert "Invalid type" in str(e.value)
    assert ".MolecularFamily" in str(e.value)
    assert ".Spectrum" in str(e.value)


def test_get_links_no_linkfinder(npl, gcfs):
    """Test `get_links` method when no LinkFinder object is found."""
    mc = MetcalfScoring(npl)
    mc.LINKFINDER = None
    with pytest.raises(ValueError) as e:
        mc.get_links(*gcfs, link_collection=LinkCollection())
    assert "LinkFinder object not found." in str(e.value)
