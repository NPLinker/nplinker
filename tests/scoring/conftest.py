from pytest import fixture
from nplinker.genomics import GCF
from nplinker.metabolomics.molecular_family import MolecularFamily
from nplinker.metabolomics.spectrum import Spectrum
from nplinker.strain_collection import StrainCollection
from nplinker.strains import Strain


@fixture(scope='session')
def strains_list() -> tuple[Strain, Strain, Strain]:
    return Strain('strain1'), Strain('strain2'), Strain('strain3')


@fixture(scope='session')
def strains(strains_list) -> StrainCollection:
    strains = StrainCollection()
    for strain in strains_list:
        strains.add(strain)
    return strains


@fixture(scope='session')
def gcfs(strains_list) -> tuple[GCF, GCF, GCF]:
    gcf1 = GCF('gcf1')
    gcf1.strains.add(strains_list[0])
    gcf2 = GCF('gcf2')
    gcf2.strains.add(strains_list[1])
    gcf3 = GCF('gcf3')
    gcf3.strains.add(strains_list[0])
    gcf3.strains.add(strains_list[1])
    return gcf1, gcf2, gcf3


@fixture(scope='session')
def spectra(strains_list) -> tuple[Spectrum, Spectrum, Spectrum]:
    spectrum1 = Spectrum(1, [(1, 1)], "spectrum1", None)
    spectrum1.strains.add(strains_list[0])
    spectrum2 = Spectrum(2, [(1, 1)], "spectrum2", None)
    spectrum2.strains.add(strains_list[1])
    spectrum3 = Spectrum(3, [(1, 1)], "spectrum3", None)
    spectrum3.strains.add(strains_list[0])
    spectrum3.strains.add(strains_list[1])
    return spectrum1, spectrum2, spectrum3


@fixture(scope='session')
def mfs(spectra) -> tuple[MolecularFamily, MolecularFamily, MolecularFamily]:
    """For simplicity, we just use one Spectrum object for each MolecularFamily
    object, and notice that they are not SingletonFamily object.
    """
    mf1 = MolecularFamily('mf1')
    mf1.add_spectrum(spectra[0])
    mf2 = MolecularFamily('mf2')
    mf2.add_spectrum(spectra[1])
    mf3 = MolecularFamily('mf3')
    mf3.add_spectrum(spectra[2])
    return mf1, mf2, mf3
