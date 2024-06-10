import pytest
from pytest import fixture
from nplinker.scoring import LinkGraph
from nplinker.scoring import Score


@fixture(scope="module")
def score():
    return Score("metcalf", 1.0, {"cutoff": 0.5})


@fixture
def lg(gcfs, spectra, score):
    lg = LinkGraph()
    lg.add_link(gcfs[0], spectra[0], metcalf=score)
    return lg


def test_init():
    lg = LinkGraph()
    assert len(lg) == 0


def test_len(lg):
    assert len(lg) == 2  # 2 objects or nodes


def test_getitem(lg, gcfs, spectra, score):
    # test existing objects
    assert lg[gcfs[0]] == {spectra[0]: {"metcalf": score}}
    assert lg[spectra[0]] == {gcfs[0]: {"metcalf": score}}

    # test non-existing object
    with pytest.raises(KeyError, match=".* not found in the link graph."):
        lg[gcfs[1]]

    # test invalid object
    with pytest.raises(TypeError, match=".* is not a GCF, Spectrum, or MolecularFamily object."):
        lg["gcf"]


def test_links(lg, gcfs, spectra, score):
    assert len(lg.links) == 1
    assert lg.links == [(gcfs[0], spectra[0], {"metcalf": score})]


def test_add_link(gcfs, spectra, score):
    lg = LinkGraph()
    lg.add_link(gcfs[0], spectra[0], metcalf=score)

    # test invalid objects
    with pytest.raises(TypeError, match=".* is not a GCF, Spectrum, or MolecularFamily object."):
        lg.add_link("gcf", spectra[0], metcalf=score)

    with pytest.raises(TypeError, match=".* is not a Spectrum or MolecularFamily object."):
        lg.add_link(gcfs[0], "spectrum", metcalf=score)

    with pytest.raises(TypeError, match=".* is not a Spectrum or MolecularFamily object."):
        lg.add_link(gcfs[0], gcfs[0], metcalf=score)

    with pytest.raises(TypeError, match=".* is not a GCF object."):
        lg.add_link(spectra[0], "gcf", metcalf=score)

    # test invalid scoring data
    with pytest.raises(
        ValueError, match="At least one scoring method and its data must be provided."
    ):
        lg.add_link(gcfs[0], spectra[0])

    with pytest.raises(ValueError, match=".* is not a valid name of scoring method.*"):
        lg.add_link(gcfs[0], spectra[0], invalid=score)

    with pytest.raises(TypeError, match=".* is not a Score object."):
        lg.add_link(gcfs[0], spectra[0], metcalf="score")


def test_has_link(lg, gcfs, spectra):
    assert lg.has_link(gcfs[0], spectra[0]) is True
    assert lg.has_link(gcfs[0], spectra[1]) is False
    assert lg.has_link(gcfs[1], spectra[1]) is False


def test_get_link_data(lg, gcfs, spectra, score):
    assert lg.get_link_data(gcfs[0], spectra[0]) == {"metcalf": score}
    assert lg.get_link_data(gcfs[0], spectra[1]) is None
