import pytest
from nplinker.metabolomics.gnps.gnps_annotation_loader import \
    GNPSAnnotationLoader
from nplinker.metabolomics.gnps.gnps_format import GNPSFormat


@pytest.mark.parametrize("workflow, expected",
                         [[GNPSFormat.FBMN, 43], [GNPSFormat.SNETS, 1154],
                          [GNPSFormat.SNETSV2, 18]])
def test_annotation_loader(workflow, expected, gnps_annotations_files):
    loader = GNPSAnnotationLoader(gnps_annotations_files[workflow])
    assert len(loader.annotations) == expected  # number of spectra

    required_columns = [
        "#Scan#", "Compound_Name", 'Organism', 'MQScore', 'SpectrumID',
        "png_url", "json_url", "svg_url", "spectrum_url"
    ]
    for annotations in loader.annotations.values():
        assert isinstance(annotations, dict)
        for name in required_columns:
            assert name in annotations
