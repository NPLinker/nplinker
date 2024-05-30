import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal
from nplinker.scoring import LinkCollection
from nplinker.scoring import MetcalfScoring
from nplinker.scoring import ObjectLink
from nplinker.scoring.linking import DataLinks


def test_init(npl):
    mc = MetcalfScoring(npl)
    assert mc.npl == npl
    assert mc.name == "metcalf"
    assert mc.cutoff == 1.0
    assert mc.standardised is True
    assert mc.DATALINKS is None
    assert_frame_equal(mc.raw_score_spec_gcf, pd.DataFrame())
    assert_frame_equal(mc.raw_score_mf_gcf, pd.DataFrame())
    assert mc.metcalf_mean is None
    assert mc.metcalf_std is None


#
# Test the `setup` method
#


def test_setup(mc, datalinks):
    """Test `setup` method when cache file does not exist."""
    assert isinstance(mc.DATALINKS, DataLinks)

    assert_frame_equal(mc.DATALINKS.occurrence_gcf_strain, datalinks.occurrence_gcf_strain)
    assert_frame_equal(mc.DATALINKS.cooccurrence_spec_gcf, datalinks.cooccurrence_spec_gcf)

    assert_frame_equal(
        mc.raw_score_spec_gcf,
        pd.DataFrame(
            [[12, -9, 11], [-9, 12, 11], [1, 1, 21]],
            index=["spectrum1", "spectrum2", "spectrum3"],
            columns=["gcf1", "gcf2", "gcf3"],
        ),
    )
    assert_frame_equal(
        mc.raw_score_mf_gcf,
        pd.DataFrame(
            [[12, -9, 11], [-9, 12, 11], [1, 1, 21]],
            index=["mf1", "mf2", "mf3"],
            columns=["gcf1", "gcf2", "gcf3"],
        ),
    )

    assert isinstance(mc.metcalf_mean, np.ndarray)
    assert isinstance(mc.metcalf_std, np.ndarray)
    assert mc.metcalf_mean.shape == (4, 4)  # (n_strains+1 , n_strains+1)
    assert mc.metcalf_mean.shape == (4, 4)


def test_setup_load_cache(mc, npl, datalinks, caplog):
    """Test `setup` method when cache file exists."""
    mc.setup(npl)

    assert isinstance(mc.DATALINKS, DataLinks)

    assert_frame_equal(mc.DATALINKS.occurrence_gcf_strain, datalinks.occurrence_gcf_strain)
    assert_frame_equal(mc.DATALINKS.cooccurrence_spec_gcf, datalinks.cooccurrence_spec_gcf)

    assert_frame_equal(
        mc.raw_score_spec_gcf,
        pd.DataFrame(
            [[12, -9, 11], [-9, 12, 11], [1, 1, 21]],
            index=["spectrum1", "spectrum2", "spectrum3"],
            columns=["gcf1", "gcf2", "gcf3"],
        ),
    )
    assert_frame_equal(
        mc.raw_score_mf_gcf,
        pd.DataFrame(
            [[12, -9, 11], [-9, 12, 11], [1, 1, 21]],
            index=["mf1", "mf2", "mf3"],
            columns=["gcf1", "gcf2", "gcf3"],
        ),
    )

    assert isinstance(mc.metcalf_mean, np.ndarray)
    assert isinstance(mc.metcalf_std, np.ndarray)
    assert mc.metcalf_mean.shape == (4, 4)  # (n_strains+1 , n_strains+1)
    assert mc.metcalf_mean.shape == (4, 4)


#
# Test the `calc_score` method
#


def test_calc_score_raw_score(mc, datalinks):
    """Test `calc_score` method for `raw_score_spec_gcf` and `raw_score_mf_gcf`.

    The expected values are calculated manually by using values from `test_init`
    of `test_data_links.py` and the default scoring weights.
    """
    # link type = 'spec-gcf'
    mc.calc_score(datalinks, link_type="spec-gcf")
    assert_frame_equal(
        mc.raw_score_spec_gcf,
        pd.DataFrame(
            [[12, -9, 11], [-9, 12, 11], [1, 1, 21]],
            index=["spectrum1", "spectrum2", "spectrum3"],
            columns=["gcf1", "gcf2", "gcf3"],
        ),
    )
    # link type = 'mf-gcf'
    mc.calc_score(datalinks, link_type="mf-gcf")
    assert_frame_equal(
        mc.raw_score_mf_gcf,
        pd.DataFrame(
            [[12, -9, 11], [-9, 12, 11], [1, 1, 21]],
            index=["mf1", "mf2", "mf3"],
            columns=["gcf1", "gcf2", "gcf3"],
        ),
    )


def test_calc_score_mean_std(mc, datalinks):
    """Test `calc_score` method for `metcalf_mean` and `metcalf_std`."""
    mc.calc_score(datalinks, link_type="spec-gcf")
    assert isinstance(mc.metcalf_mean, np.ndarray)
    assert isinstance(mc.metcalf_std, np.ndarray)
    assert mc.metcalf_mean.shape == (4, 4)  # (n_strains+1 , n_strains+1)
    assert mc.metcalf_mean.shape == (4, 4)
    # TODO CG: add tests for values after refactoring _calc_mean_std method
    # assert mc.metcalf_mean == expected_array


#
# Test the `get_links` method
#


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
    assert {i.gcf_id for i in links.keys()} == {"gcf1", "gcf2", "gcf3"}
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
    assert {i.gcf_id for i in links.keys()} == {"gcf1", "gcf2", "gcf3"}
    assert isinstance(links[gcfs[0]][spectra[0]], ObjectLink)
    assert links[gcfs[0]][spectra[0]].data(mc) == 12
    assert links[gcfs[1]].get(spectra[0]) is None
    assert links[gcfs[2]][spectra[0]].data(mc) == 11
    assert links[gcfs[0]][mfs[0]].data(mc) == 12
    assert links[gcfs[1]][mfs[1]].data(mc) == 12
    assert links[gcfs[2]][mfs[2]].data(mc) == 21


@pytest.mark.skip(reason="To add after refactoring relevant code.")
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
    assert {i.spectrum_id for i in links.keys()} == {"spectrum1", "spectrum2", "spectrum3"}
    assert isinstance(links[spectra[0]][gcfs[0]], ObjectLink)
    assert links[spectra[0]][gcfs[0]].data(mc) == 12
    assert links[spectra[0]][gcfs[1]].data(mc) == -9
    assert links[spectra[0]][gcfs[2]].data(mc) == 11

    mc.cutoff = 0
    links = mc.get_links(*spectra, link_collection=LinkCollection())
    assert isinstance(links, LinkCollection)
    links = links.links  # dict of link values
    assert len(links) == 3
    assert {i.spectrum_id for i in links.keys()} == {"spectrum1", "spectrum2", "spectrum3"}
    assert isinstance(links[spectra[0]][gcfs[0]], ObjectLink)
    assert links[spectra[0]][gcfs[0]].data(mc) == 12
    assert links[spectra[0]].get(gcfs[1]) is None
    assert links[spectra[0]][gcfs[2]].data(mc) == 11


@pytest.mark.skip(reason="To add after refactoring relevant code.")
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
    assert {i.family_id for i in links.keys()} == {"mf1", "mf2", "mf3"}
    assert isinstance(links[mfs[0]][gcfs[0]], ObjectLink)
    assert links[mfs[0]][gcfs[0]].data(mc) == 12
    assert links[mfs[0]][gcfs[1]].data(mc) == -9
    assert links[mfs[0]][gcfs[2]].data(mc) == 11

    mc.cutoff = 0
    links = mc.get_links(*mfs, link_collection=LinkCollection())
    assert isinstance(links, LinkCollection)
    links = links.links
    assert len(links) == 3
    assert {i.family_id for i in links.keys()} == {"mf1", "mf2", "mf3"}
    assert isinstance(links[mfs[0]][gcfs[0]], ObjectLink)
    assert links[mfs[0]][gcfs[0]].data(mc) == 12
    assert links[mfs[0]].get(gcfs[1]) is None
    assert links[mfs[0]][gcfs[2]].data(mc) == 11


@pytest.mark.skip(reason="To add after refactoring relevant code.")
def test_get_links_mf_standardised_true(mc, gcfs, mfs):
    """Test `get_links` method when input is MolecularFamily objects and `standardised` is True."""
    mc.standardised = True
    ...


@pytest.mark.parametrize(
    "objects, expected", [([], "Empty input objects"), ("", "Empty input objects")]
)
def test_get_links_invalid_input_value(mc, objects, expected):
    with pytest.raises(ValueError) as e:
        mc.get_links(*objects, link_collection=LinkCollection())
    assert expected in str(e.value)


@pytest.mark.parametrize(
    "objects, expected",
    [
        ([1], "Invalid type {<class 'int'>}"),
        ([1, 2], "Invalid type {<class 'int'>}"),
        ("12", "Invalid type {<class 'str'>}"),
    ],
)
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


#
# Test the `_get_links` method
#


def test__get_links_gcf(mc, datalinks, gcfs):
    """Test `get_links` method for input GCF objects."""
    mc.calc_score(datalinks, link_type="spec-gcf")
    mc.calc_score(datalinks, link_type="mf-gcf")
    index_names = ["source", "target", "score"]

    # cutoff = negative infinity (float)
    links = mc._get_links(*gcfs, score_cutoff=np.NINF)
    assert len(links) == 2
    # expected values got from `test_calc_score_raw_score`
    assert_frame_equal(
        links[0],
        pd.DataFrame(
            [
                ["gcf1", "gcf2", "gcf3"] * 3,
                [
                    *["spectrum1"] * 3,
                    *["spectrum2"] * 3,
                    *["spectrum3"] * 3,
                ],
                [12, -9, 11, -9, 12, 11, 1, 1, 21],
            ],
            index=index_names,
        ),
    )
    assert_frame_equal(
        links[1],
        pd.DataFrame(
            [
                ["gcf1", "gcf2", "gcf3"] * 3,
                [
                    *["mf1"] * 3,
                    *["mf2"] * 3,
                    *["mf3"] * 3,
                ],
                [12, -9, 11, -9, 12, 11, 1, 1, 21],
            ],
            index=index_names,
        ),
    )

    # cutoff = 0
    links = mc._get_links(*gcfs, score_cutoff=0)
    assert len(links) == 2
    assert_frame_equal(
        links[0],
        pd.DataFrame(
            [
                ["gcf1", "gcf3", "gcf2", "gcf3", "gcf1", "gcf2", "gcf3"],
                [
                    *["spectrum1"] * 2,
                    *["spectrum2"] * 2,
                    *["spectrum3"] * 3,
                ],
                [12, 11, 12, 11, 1, 1, 21],
            ],
            index=index_names,
        ),
    )
    assert_frame_equal(
        links[1],
        pd.DataFrame(
            [
                ["gcf1", "gcf3", "gcf2", "gcf3", "gcf1", "gcf2", "gcf3"],
                [
                    *["mf1"] * 2,
                    *["mf2"] * 2,
                    *["mf3"] * 3,
                ],
                [12, 11, 12, 11, 1, 1, 21],
            ],
            index=index_names,
        ),
    )


def test__get_links_spec(mc, datalinks, spectra):
    """Test `get_links` method for input Spectrum objects."""
    mc.calc_score(datalinks, link_type="spec-gcf")
    mc.calc_score(datalinks, link_type="mf-gcf")
    index_names = ["source", "target", "score"]
    # cutoff = negative infinity (float)
    links = mc._get_links(*spectra, score_cutoff=np.NINF)
    assert len(links) == 1
    assert_frame_equal(
        links[0],
        pd.DataFrame(
            [
                [
                    *["spectrum1"] * 3,
                    *["spectrum2"] * 3,
                    *["spectrum3"] * 3,
                ],
                ["gcf1", "gcf2", "gcf3"] * 3,
                [12, -9, 11, -9, 12, 11, 1, 1, 21],
            ],
            index=index_names,
        ),
    )
    # cutoff = 0
    links = mc._get_links(*spectra, score_cutoff=0)
    assert_frame_equal(
        links[0],
        pd.DataFrame(
            [
                [
                    *["spectrum1"] * 2,
                    *["spectrum2"] * 2,
                    *["spectrum3"] * 3,
                ],
                ["gcf1", "gcf3", "gcf2", "gcf3", "gcf1", "gcf2", "gcf3"],
                [12, 11, 12, 11, 1, 1, 21],
            ],
            index=index_names,
        ),
    )


def test__get_links_mf(mc, datalinks, mfs):
    """Test `get_links` method for input MolecularFamily objects."""
    mc.calc_score(datalinks, link_type="spec-gcf")
    mc.calc_score(datalinks, link_type="mf-gcf")
    index_names = ["source", "target", "score"]
    # cutoff = negative infinity (float)
    links = mc._get_links(*mfs, score_cutoff=np.NINF)
    assert len(links) == 1
    assert_frame_equal(
        links[0],
        pd.DataFrame(
            [
                [
                    *["mf1"] * 3,
                    *["mf2"] * 3,
                    *["mf3"] * 3,
                ],
                ["gcf1", "gcf2", "gcf3"] * 3,
                [12, -9, 11, -9, 12, 11, 1, 1, 21],
            ],
            index=index_names,
        ),
    )
    # cutoff = 0
    links = mc._get_links(*mfs, score_cutoff=0)
    assert_frame_equal(
        links[0],
        pd.DataFrame(
            [
                [
                    *["mf1"] * 2,
                    *["mf2"] * 2,
                    *["mf3"] * 3,
                ],
                ["gcf1", "gcf3", "gcf2", "gcf3", "gcf1", "gcf2", "gcf3"],
                [12, 11, 12, 11, 1, 1, 21],
            ],
            index=index_names,
        ),
    )


@pytest.mark.parametrize(
    "objects, expected", [([], "Empty input objects"), ("", "Empty input objects")]
)
def test_get_links_invalid_value(mc, objects, expected):
    with pytest.raises(ValueError) as e:
        mc._get_links(*objects)
    assert expected in str(e.value)


@pytest.mark.parametrize(
    "objects, expected",
    [
        ([1], "Invalid type {<class 'int'>}"),
        ([1, 2], "Invalid type {<class 'int'>}"),
        ("12", "Invalid type {<class 'str'>}"),
    ],
)
def test__get_links_invalid_type(mc, objects, expected):
    with pytest.raises(TypeError) as e:
        mc._get_links(*objects)
    assert expected in str(e.value)


def test__get_links_invalid_mixed_types(mc, spectra, mfs):
    objects = (*spectra, *mfs)
    with pytest.raises(TypeError) as e:
        mc._get_links(*objects)
    assert "Invalid type" in str(e.value)
    assert ".MolecularFamily" in str(e.value)
    assert ".Spectrum" in str(e.value)
