import os

from nplinker.annotations import load_annotations
from .test_metabolomics import spec_dict

testdata_dir = os.path.join(os.getcwd(),"tests", "data")



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