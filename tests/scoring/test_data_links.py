import pandas as pd
from pandas.util.testing import assert_frame_equal
from pytest import fixture
from nplinker.genomics import GCF
from nplinker.metabolomics.molecular_family import MolecularFamily
from nplinker.metabolomics.spectrum import Spectrum
from nplinker.scoring.linking import DataLinks
from nplinker.strain_collection import StrainCollection
from nplinker.strains import Strain


@fixture(scope='module')
def strains_list():
    return Strain('strain1'), Strain('strain2'), Strain('strain3')


@fixture(scope='module')
def strains(strains_list):
    strains = StrainCollection()
    for strain in strains_list:
        strains.add(strain)
    return strains


@fixture(scope='module')
def gcfs(strains_list):
    gcf1 = GCF('gcf1')
    gcf1.strains.add(strains_list[0])
    gcf2 = GCF('gcf2')
    gcf2.strains.add(strains_list[1])
    gcf3 = GCF('gcf3')
    gcf3.strains.add(strains_list[0])
    gcf3.strains.add(strains_list[1])
    return gcf1, gcf2, gcf3


@fixture(scope='module')
def spectra(strains_list):
    spectrum1 = Spectrum(1, [(1, 1)], "spectrum1", None)
    spectrum1.strains.add(strains_list[0])
    spectrum2 = Spectrum(2, [(1, 1)], "spectrum2", None)
    spectrum2.strains.add(strains_list[1])
    spectrum3 = Spectrum(3, [(1, 1)], "spectrum3", None)
    spectrum3.strains.add(strains_list[0])
    spectrum3.strains.add(strains_list[1])
    return spectrum1, spectrum2, spectrum3


@fixture(scope='module')
def mfs(spectra):
    mf1 = MolecularFamily('mf1')
    mf1.add_spectrum(spectra[0])
    mf2 = MolecularFamily('mf2')
    mf2.add_spectrum(spectra[1])
    mf3 = MolecularFamily('mf3')
    mf3.add_spectrum(spectra[2])
    return mf1, mf2, mf3


@fixture(scope='module')
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


# TODO: add tests for the DataLinks class
