

import os
from nplinker.metabolomics.load_gnps import GNPS_FORMAT_NEW_FBMN, GNPS_FORMAT_OLD_ALLFILES, _identify_gnps_format, _load_clusterinfo_old, load_gnps
from nplinker.strain_collection import StrainCollection

from .test_metabolomics import spec_dict

testdata_dir = os.path.join(os.getcwd(),"tests", "data")
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

    assert len(unknown_strains) is not None


def test_identify_gnps_format():
    actual = _identify_gnps_format(nodes_file, None)

    assert actual is not GNPS_FORMAT_NEW_FBMN
    assert actual is GNPS_FORMAT_OLD_ALLFILES


def test_load_clusterinfo_old(spec_dict):
    sut = _load_clusterinfo_old(
        GNPS_FORMAT_OLD_ALLFILES,
        strains,
        nodes_file,
        spec_dict
    )
    
    assert len(sut) is not None
