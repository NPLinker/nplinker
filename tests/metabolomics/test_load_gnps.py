from itertools import chain

import pytest
from nplinker.metabolomics.load_gnps import _messy_strain_naming_lookup, _parse_mzxml_header
from nplinker.metabolomics.gnps.gnps_format import gnps_format_from_file_mapping, GNPSFormat
from nplinker.metabolomics.load_gnps import _load_clusterinfo_old
from nplinker.metabolomics.load_gnps import load_gnps
from nplinker.strain_collection import StrainCollection
from nplinker.utils import get_headers

from .. import DATA_DIR


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


def test_load_clusterinfo_old(spec_dict):
    metadata_keys_before = set(chain(*[x.metadata.copy().keys() for x in spec_dict.values()]))

    assert 'files' not in metadata_keys_before
    assert 'cluster_index' not in metadata_keys_before


    sut = _load_clusterinfo_old(
        GNPSFormat.AllFiles,
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
    ["blub", None],
    ["42b.mzXML.copy", "42b.mzXML"],
    ["Salinispora arenicola CNB527_blub", "42b.mzXML"],
    ["GCF_000514775.1", "9b.mzXML"]
])
def test_messy_strain_naming_lookup(collection_from_file: StrainCollection, messy_alias: str, expected: str|None):
    actual = _messy_strain_naming_lookup(messy_alias, collection_from_file)

    if expected is not None:
        assert actual == collection_from_file.lookup(expected)
    else:
        assert actual == expected


def test_parse_mzxml_header():
    headers = get_headers(str(DATA_DIR / "nodes_fbmn.csv"))
    hdr = headers[10]
    actual = _parse_mzxml_header(hdr, StrainCollection(), None, None)
    assert actual is not None
