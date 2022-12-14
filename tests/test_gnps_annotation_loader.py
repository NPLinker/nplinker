import pytest
from typing_extensions import Self
from nplinker.annotations import load_annotations
from nplinker.metabolomics.gnps.gnps_annotation_loader import \
    GNPSAnnotationLoader
from nplinker.metabolomics.spectrum import Spectrum
from . import DATA_DIR
from .test_metabolomics import spec_dict


class GNPSAnnotationLoaderBuilder:
    def __init__(self):
        self._file = DATA_DIR / "ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-92036537-download_cytoscape_data/DB_result/7dc5b46b50d94246a1de12ef485d0f75.tsv"
    
    def with_file(self, file) -> Self:
        self._file = file
        return self
    
    def build(self) -> GNPSAnnotationLoader:
        return GNPSAnnotationLoader(self._file)

def test_default():
    sut = GNPSAnnotationLoaderBuilder().build()
    assert sut is not None


@pytest.mark.parametrize("file, expected", [
    [DATA_DIR / "ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-92036537-download_cytoscape_data/DB_result/7dc5b46b50d94246a1de12ef485d0f75.tsv", 43],
    [DATA_DIR / "ProteoSAFe-METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra/result_specnets_DB/885e4c5485ba42569e4876d1fe90d759.tsv", 1154]
])
def test_reads_all_annotations(file, expected):
    sut = GNPSAnnotationLoaderBuilder().with_file(file).build()
    assert len(sut.annotations()) == expected


def test_annotations_are_equal(spec_dict: dict[int, Spectrum]):
    annotations_dir = DATA_DIR / "annotations"
    annotations_file = annotations_dir / "gnps_annotations.tsv"

    spectra = list(spec_dict.values())

    sut = load_annotations(
        annotations_dir,
        "",
        spectra,
        spec_dict
    )
    expected: dict[int, dict] = {}
    for x in sut:
        if x.has_annotations():
            expected[x.spectrum_id] = x.gnps_annotations

    actual = GNPSAnnotationLoaderBuilder().with_file(annotations_file).build().annotations()
    
    for key in expected.keys():
        assert expected[key].items() <= actual[key].items()