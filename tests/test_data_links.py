from nplinker.scoring.linking.data_linking import DataLinks
from .test_metabolomics import spec_dict, spec_with_families
from .test_gnps_molecular_family_loader import molecular_families_gnps


def test_collect_mappings_spec(spec_with_families):
    sut = DataLinks()
    sut.collect_mappings_spec(spec_with_families.values())
    actual = sut.mapping_spec.shape

    assert actual == (25935,3)


def test_collect_mappings_spec_v2(molecular_families_gnps):
    sut = DataLinks()
    sut.collect_mappings_spec_v2(molecular_families_gnps)
    actual = sut.mapping_spec.shape

    assert actual == (25935,3)


def test_mappings_are_equal(spec_with_families, molecular_families_gnps):
    sut = DataLinks()
    sut.collect_mappings_spec(spec_with_families.values())
    actual = sut.mapping_spec

    sut.collect_mappings_spec_v2(molecular_families_gnps)
    expected = sut.mapping_spec

    assert actual.eq(expected).all(axis=None)