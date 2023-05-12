import numpy as np
import pytest
from pytest import fixture
from nplinker.nplinker import NPLinker
from nplinker.scoring import LinkCollection
from nplinker.scoring import MetcalfScoring
from nplinker.scoring import ObjectLink
from nplinker.scoring.linking import DataLinks
from nplinker.scoring.linking import LinkFinder
from .. import DATA_DIR


@fixture(scope='module')
def datalinks(gcfs, spectra, mfs, strains) -> DataLinks:
    """DataLinks object. See `test_data_links.py` for its values."""
    return DataLinks(gcfs, spectra, mfs, strains)


@fixture(scope='module')
def linkfinder(datalinks) -> LinkFinder:
    """LinkFinder object. See `test_link_finder.py` for its values."""
    linkfinder = LinkFinder()
    linkfinder.cal_score(datalinks, link_type='spec-gcf')
    linkfinder.cal_score(datalinks, link_type='mf-gcf')
    return linkfinder


@fixture(scope='module')
def npl(gcfs, spectra, mfs, strains, tmp_path_factory) -> NPLinker:
    """Constructed NPLinker object.

    This NPLinker object does not do loading `npl.load_data()`, instead we
    manually set its attributes to the values we want to test.

    The config file `nplinker_demo1.toml` does not affect the tests, just
    making sure the NPLinker object can be created succesfully.
    """
    npl = NPLinker(str(DATA_DIR / 'nplinker_demo1.toml'))
    npl._gcfs = gcfs
    npl._spectra = spectra
    npl._molfams = mfs
    npl._strains = strains
    npl._gcf_lookup = {gcf.gcf_id: gcf for gcf in gcfs}
    npl._mf_lookup = {mf.family_id: mf for mf in mfs}
    npl._spec_lookup = {spec.spectrum_id: spec for spec in spectra}
    # tmp path to store 'metcalf/metcalf_scores.pckl' file
    # Must use `tmp_path_factory` (session scope) instead of `tmp_path` (function scope)
    npl._loader._root = tmp_path_factory.mktemp('npl_test')
    return npl


@fixture(scope='module')
def mc(npl) -> MetcalfScoring:
    """MetcalfScoring object."""
    mc = MetcalfScoring(npl)
    mc.setup(npl)
    return mc


def test_get_links_gcf_standardised_false(npl, mc, gcfs, spectra, mfs,
                                          strains_list):
    """Test `get_links` method when input is GCF objects and `standardised` is False."""
    # test raw scores (no standardisation)
    mc.standardised = False

    # when cutoff is negative infinity, i.e. taking all scores
    mc.cutoff = np.NINF
    links = npl.get_links(list(gcfs), mc, and_mode=True)
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
    # expected values are from `test_get_common_strains_spec` of test_data_links.py
    assert links[gcfs[0]][spectra[0]].shared_strains == [strains_list[0]]
    assert links[gcfs[1]][spectra[0]].shared_strains == []
    assert links[gcfs[2]][spectra[0]].shared_strains == [strains_list[0]]
    assert links[gcfs[0]][mfs[0]].shared_strains == [strains_list[0]]
    assert links[gcfs[1]][mfs[1]].shared_strains == [strains_list[1]]
    assert set(links[gcfs[2]][mfs[2]].shared_strains) == set(strains_list[0:2])

    # when test cutoff is 0, i.e. taking scores >= 0
    mc.cutoff = 0
    links = npl.get_links(list(gcfs), mc, and_mode=True)
    assert isinstance(links, LinkCollection)
    links = links.links
    assert {i.gcf_id for i in links.keys()} == {'gcf1', 'gcf2', 'gcf3'}
    assert isinstance(links[gcfs[0]][spectra[0]], ObjectLink)
    # test scores
    assert links[gcfs[0]][spectra[0]].data(mc) == 12
    assert links[gcfs[1]].get(spectra[0]) is None
    assert links[gcfs[2]][spectra[0]].data(mc) == 11
    assert links[gcfs[0]][mfs[0]].data(mc) == 12
    assert links[gcfs[1]][mfs[1]].data(mc) == 12
    assert links[gcfs[2]][mfs[2]].data(mc) == 21
    # test shared strains
    assert links[gcfs[0]][spectra[0]].shared_strains == [strains_list[0]]
    assert links[gcfs[2]][spectra[0]].shared_strains == [strains_list[0]]
    assert links[gcfs[0]][mfs[0]].shared_strains == [strains_list[0]]
    assert links[gcfs[1]][mfs[1]].shared_strains == [strains_list[1]]
    assert set(links[gcfs[2]][mfs[2]].shared_strains) == set(strains_list[0:2])


@pytest.mark.skip(reason='To add after refactoring relevant code.')
def test_get_links_gcf_standardised_true(npl, mc, gcfs, spectra, mfs,
                                         strains_list):
    """Test `get_links` method when input is GCF objects and `standardised` is True."""
    mc.standardised = True
    ...


def test_get_links_spec_standardised_false(npl, mc, gcfs, spectra,
                                           strains_list):
    """Test `get_links` method when input is Spectrum objects and `standardised` is False."""
    mc.standardised = False

    mc.cutoff = np.NINF
    links = npl.get_links(list(spectra), mc, and_mode=True)
    assert isinstance(links, LinkCollection)
    links = links.links  # dict of link values
    assert len(links) == 3
    assert {i.spectrum_id
            for i in links.keys()} == {'spectrum1', 'spectrum2', 'spectrum3'}
    assert isinstance(links[spectra[0]][gcfs[0]], ObjectLink)
    assert links[spectra[0]][gcfs[0]].data(mc) == 12
    assert links[spectra[0]][gcfs[1]].data(mc) == -9
    assert links[spectra[0]][gcfs[2]].data(mc) == 11
    assert links[spectra[0]][gcfs[0]].shared_strains == [strains_list[0]]
    assert links[spectra[0]][gcfs[1]].shared_strains == []
    assert links[spectra[0]][gcfs[2]].shared_strains == [strains_list[0]]

    mc.cutoff = 0
    links = npl.get_links(list(spectra), mc, and_mode=True)
    assert isinstance(links, LinkCollection)
    links = links.links  # dict of link values
    assert len(links) == 3
    assert {i.spectrum_id
            for i in links.keys()} == {'spectrum1', 'spectrum2', 'spectrum3'}
    assert isinstance(links[spectra[0]][gcfs[0]], ObjectLink)
    assert links[spectra[0]][gcfs[0]].data(mc) == 12
    assert links[spectra[0]].get(gcfs[1]) is None
    assert links[spectra[0]][gcfs[2]].data(mc) == 11
    assert links[spectra[0]][gcfs[0]].shared_strains == [strains_list[0]]
    assert links[spectra[0]][gcfs[2]].shared_strains == [strains_list[0]]


@pytest.mark.skip(reason='To add after refactoring relevant code.')
def test_get_links_spec_standardised_true(npl, mc, gcfs, spectra,
                                          strains_list):
    """Test `get_links` method when input is Spectrum objects and `standardised` is True."""
    mc.standardised = True
    ...


def test_get_links_mf_standardised_false(npl, mc, gcfs, mfs, strains_list):
    """Test `get_links` method when input is MolecularFamily objects and `standardised` is False."""
    mc.standardised = False

    mc.cutoff = np.NINF
    links = npl.get_links(list(mfs), mc, and_mode=True)
    assert isinstance(links, LinkCollection)
    links = links.links
    assert len(links) == 3
    assert {i.family_id for i in links.keys()} == {'mf1', 'mf2', 'mf3'}
    assert isinstance(links[mfs[0]][gcfs[0]], ObjectLink)
    assert links[mfs[0]][gcfs[0]].data(mc) == 12
    assert links[mfs[0]][gcfs[1]].data(mc) == -9
    assert links[mfs[0]][gcfs[2]].data(mc) == 11
    assert links[mfs[0]][gcfs[0]].shared_strains == [strains_list[0]]
    assert links[mfs[0]][gcfs[1]].shared_strains == []
    assert links[mfs[0]][gcfs[2]].shared_strains == [strains_list[0]]

    mc.cutoff = 0
    links = npl.get_links(list(mfs), mc, and_mode=True)
    assert isinstance(links, LinkCollection)
    links = links.links
    assert len(links) == 3
    assert {i.family_id for i in links.keys()} == {'mf1', 'mf2', 'mf3'}
    assert isinstance(links[mfs[0]][gcfs[0]], ObjectLink)
    assert links[mfs[0]][gcfs[0]].data(mc) == 12
    assert links[mfs[0]].get(gcfs[1]) is None
    assert links[mfs[0]][gcfs[2]].data(mc) == 11
    assert links[mfs[0]][gcfs[0]].shared_strains == [strains_list[0]]
    assert links[mfs[0]][gcfs[2]].shared_strains == [strains_list[0]]


@pytest.mark.skip(reason='To add after refactoring relevant code.')
def test_get_links_mf_standardised_true(npl, mc, gcfs, mfs, strains_list):
    """Test `get_links` method when input is MolecularFamily objects and `standardised` is True."""
    mc.standardised = True
    ...
