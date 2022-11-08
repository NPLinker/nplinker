from itertools import chain
from typing import Optional

import pytest
from nplinker.metabolomics.load_gnps import GNPS_FORMAT_NEW_FBMN, _messy_strain_naming_lookup
from nplinker.metabolomics.load_gnps import GNPS_FORMAT_OLD_ALLFILES
from nplinker.metabolomics.load_gnps import identify_gnps_format
from nplinker.metabolomics.load_gnps import _load_clusterinfo_old
from nplinker.metabolomics.load_gnps import load_gnps
from nplinker.strain_collection import StrainCollection
from .test_metabolomics import spec_dict
from .test_strain_collection import collection_from_file

from . import DATA_DIR


nodes_file = DATA_DIR / "nodes.tsv"
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


@pytest.mark.parametrize("filename, gnps_format", [
    [nodes_file, GNPS_FORMAT_OLD_ALLFILES],
    [ DATA_DIR / "nodes_fbmn.csv", GNPS_FORMAT_NEW_FBMN]
])
def test_identify_gnps_format(filename, gnps_format):
    actual = identify_gnps_format(filename, None)

    assert actual is gnps_format


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


@pytest.mark.parametrize("messy_alias, expected", [
    ["42b.mzXML", "42b.mzXML"],
    ["42b.mzXML.copy", "42b.mzXML"],
    ["blub", None],
    ["Salinispora arenicola CNB527_blub", "42b.mzXML"],
    ["CNB527", None],
    ["GCF_000514775.1", "9b.mzXML"]
])
def test_messy_strain_naming_lookup(collection_from_file: StrainCollection, messy_alias: str, expected: Optional[str]):
    actual = _messy_strain_naming_lookup(messy_alias, collection_from_file)

    assert actual == collection_from_file.lookup(expected)
