import os
from itertools import chain
from nplinker.metabolomics.load_gnps import GNPS_FORMAT_NEW_FBMN
from nplinker.metabolomics.load_gnps import GNPS_FORMAT_OLD_ALLFILES
from nplinker.metabolomics.load_gnps import _identify_gnps_format
from nplinker.metabolomics.load_gnps import _load_clusterinfo_old
from nplinker.metabolomics.load_gnps import load_gnps
from nplinker.strain_collection import StrainCollection
from .test_metabolomics import spec_dict


testdata_dir = os.path.join(os.getcwd(), "tests", "data")
nodes_file = os.path.join(testdata_dir, "nodes.tsv")
strains = StrainCollection()


def test_load_gnps(spec_dict):
    unknown_strains = load_gnps(
        strains,
        nodes_file,
        None,
        None,
        None,
        spec_dict
    )

    assert len(unknown_strains) > 0


def test_identify_gnps_format():
    actual = _identify_gnps_format(nodes_file, None)

    assert actual is not GNPS_FORMAT_NEW_FBMN
    assert actual is GNPS_FORMAT_OLD_ALLFILES


def test_load_clusterinfo_old(spec_dict):
    metadata_keys_before = set(chain(*[x.metadata.copy().keys() for x in spec_dict.values()]))

    assert 'files' not in metadata_keys_before
    assert 'cluster_index' not in metadata_keys_before


    sut = _load_clusterinfo_old(
        GNPS_FORMAT_OLD_ALLFILES,
        strains,
        nodes_file,
        spec_dict
    )

    metadata_keys_after = set(chain(*[x.metadata.keys() for x in spec_dict.values()]))

    assert len(metadata_keys_after) == len(metadata_keys_before) + 2
  
    assert len(sut) > 0

    for spectrum_id, spec in spec_dict.items():
        metadata = spec.metadata
        assert len(metadata.get('files')) > 1
        assert isinstance(metadata.get('cluster_index'), int)
        assert spectrum_id == metadata.get('cluster_index')
