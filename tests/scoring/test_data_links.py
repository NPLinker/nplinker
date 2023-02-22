import os
import pytest
from nplinker.metabolomics.gnps.gnps_molecular_family_loader import \
    GNPSMolecularFamilyLoader
from nplinker.metabolomics.metabolomics import make_families
from nplinker.metabolomics.spectrum import Spectrum
from nplinker.scoring.linking.data_linking import DataLinks
from .. import DATA_DIR


@pytest.fixture
def spec_with_families(spec_dict) -> dict[int, Spectrum]:
    make_families(spec_dict.values())
    return spec_dict

@pytest.fixture
def molecular_families_gnps():
    filename = os.path.join(DATA_DIR, "edges.pairsinfo")
    sut = GNPSMolecularFamilyLoader(filename)
    return sut.families()

def test_collect_mappings_from_spectra(spec_with_families):
    sut = DataLinks()
    actual = sut._collect_mappings_from_spectra(spec_with_families.values())

    assert actual.shape == (25935,3)


def test_collect_mappings_from_molecular_families(molecular_families_gnps):
    sut = DataLinks()
    actual = sut._collect_mappings_from_molecular_families(molecular_families_gnps)

    assert actual.shape == (25935,3)


def test_mappings_are_equal(spec_with_families, molecular_families_gnps):
    sut = DataLinks()
    sut._collect_mappings_from_spectra(spec_with_families.values())
    actual = sut.mapping_spec

    sut._collect_mappings_from_molecular_families(molecular_families_gnps)
    expected = sut.mapping_spec

    assert actual.eq(expected).all(axis=None)
