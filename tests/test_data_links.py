from nplinker.scoring.linking.data_linking import DataLinks
from .test_metabolomics import spec_dict


def test_collect_mappings_spec(spec_dict):
    sut = DataLinks()
    sut.collect_mappings_spec(spec_dict.values())
    actual = sut.mapping_spec

    assert actual is not None

