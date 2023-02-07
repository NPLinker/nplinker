from nplinker.scoring.linking.data_linking import DataLinks
from .test_metabolomics import spec_dict, spec_with_families
from .test_gnps_molecular_family_loader import molecular_families_gnps


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