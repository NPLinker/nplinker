from pytest import fixture
from nplinker.genomics import GCF
from nplinker.metabolomics import MolecularFamily
from nplinker.metabolomics import Spectrum
from nplinker.nplinker import NPLinker
from nplinker.scoring import MetcalfScoring
from nplinker.scoring.linking import DataLinks
from nplinker.scoring.linking import LinkFinder
from nplinker.strain import Strain
from nplinker.strain_collection import StrainCollection
from .. import DATA_DIR


@fixture(scope="session")
def strains_list() -> tuple[Strain, Strain, Strain]:
    return Strain("strain1"), Strain("strain2"), Strain("strain3")


@fixture(scope="session")
def strains(strains_list) -> StrainCollection:
    strains = StrainCollection()
    for strain in strains_list:
        strains.add(strain)
    return strains


@fixture(scope="session")
def gcfs(strains_list) -> tuple[GCF, GCF, GCF]:
    gcf1 = GCF("gcf1")
    gcf1.strains.add(strains_list[0])
    gcf2 = GCF("gcf2")
    gcf2.strains.add(strains_list[1])
    gcf3 = GCF("gcf3")
    gcf3.strains.add(strains_list[0])
    gcf3.strains.add(strains_list[1])
    return gcf1, gcf2, gcf3


@fixture(scope="session")
def spectra(strains_list) -> tuple[Spectrum, Spectrum, Spectrum]:
    spectrum1 = Spectrum("spectrum1", [(1, 1)], None)
    spectrum1.strains.add(strains_list[0])
    spectrum2 = Spectrum("spectrum2", [(1, 1)], None)
    spectrum2.strains.add(strains_list[1])
    spectrum3 = Spectrum("spectrum3", [(1, 1)], None)
    spectrum3.strains.add(strains_list[0])
    spectrum3.strains.add(strains_list[1])
    return spectrum1, spectrum2, spectrum3


@fixture(scope="session")
def mfs(spectra) -> tuple[MolecularFamily, MolecularFamily, MolecularFamily]:
    """For simplicity, we just use one Spectrum object for each MolecularFamily
    object, and notice that they are not SingletonFamily object.
    """
    mf1 = MolecularFamily("mf1")
    mf1.add_spectrum(spectra[0])
    mf2 = MolecularFamily("mf2")
    mf2.add_spectrum(spectra[1])
    mf3 = MolecularFamily("mf3")
    mf3.add_spectrum(spectra[2])
    return mf1, mf2, mf3


@fixture(scope="module")
def datalinks(gcfs, spectra, mfs, strains) -> DataLinks:
    """DataLinks object. See `test_data_links.py` for its values."""
    return DataLinks(gcfs, spectra, mfs, strains)


@fixture(scope="module")
def linkfinder(datalinks) -> LinkFinder:
    """LinkFinder object. See `test_link_finder.py` for its values."""
    linkfinder = LinkFinder()
    linkfinder.calc_score(datalinks, link_type="spec-gcf")
    linkfinder.calc_score(datalinks, link_type="mf-gcf")
    return linkfinder


@fixture(scope="module")
def npl(gcfs, spectra, mfs, strains, tmp_path_factory) -> NPLinker:
    """Constructed NPLinker object.

    This NPLinker object does not do loading `npl.load_data()`, instead we
    manually set its attributes to the values we want to test.

    The config file `nplinker_demo1.toml` does not affect the tests, just
    making sure the NPLinker object can be created succesfully.
    """
    npl = NPLinker(str(DATA_DIR / "nplinker_demo1.toml"))
    npl._gcfs = gcfs
    npl._spectra = spectra
    npl._molfams = mfs
    npl._strains = strains
    npl._gcf_lookup = {gcf.gcf_id: gcf for gcf in gcfs}
    npl._mf_lookup = {mf.family_id: mf for mf in mfs}
    npl._spec_lookup = {spec.spectrum_id: spec for spec in spectra}
    # tmp path to store 'metcalf/metcalf_scores.pckl' file
    # Must use `tmp_path_factory` (session scope) instead of `tmp_path` (function scope)
    npl._loader._root = tmp_path_factory.mktemp("npl_test")
    return npl


@fixture(scope="module")
def mc(npl) -> MetcalfScoring:
    """MetcalfScoring object."""
    mc = MetcalfScoring(npl)
    mc.setup(npl)
    return mc
