import os

import pytest
from nplinker.annotations import load_annotations
from nplinker.metabolomics.metabolomics import load_spectra


testdata_dir = os.path.join(os.getcwd(),"tests", "data")


@pytest.fixture
def spec_dict() -> dict:
    mgf_file = os.path.join(testdata_dir, "spectra.mgf")
    edges_file = os.path.join(testdata_dir, "edges.pairsinfo")
    return load_spectra(mgf_file, edges_file)



def test_load_annotations(spec_dict):
    annotations_dir = os.path.join(testdata_dir,"annotations")
    annotations_file = os.path.join(annotations_dir, "annotations.tsv")

    spectra = list(spec_dict.values())

    sut = load_annotations(
        annotations_dir,
        annotations_file,
        spectra,
        spec_dict
    )

    annotations = list(filter(lambda x: len(x) > 0, [x.annotations for x in sut]))
    
    assert len(annotations) > 0
    assert len(spectra) > 0
    assert len(spec_dict.keys()) > 0