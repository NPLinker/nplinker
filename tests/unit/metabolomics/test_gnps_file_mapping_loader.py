import pytest
from nplinker.metabolomics.gnps import GNPSFileMappingLoader
from nplinker.metabolomics.gnps import GNPSFormat


@pytest.mark.parametrize(
    "workflow, num_spectra, filename",
    [
        [GNPSFormat.FBMN, 1492, "5434_5433_mod.mzXML"],
        [GNPSFormat.SNETS, 25935, "26c.mzXML"],
        [GNPSFormat.SNETSV2, 7383, "140221_ME_14_13.mzML"],
    ],
)
def test_file_mapping_loader(workflow, num_spectra, filename, gnps_file_mappings_files):
    loader = GNPSFileMappingLoader(gnps_file_mappings_files[workflow])
    assert len(loader.mappings) == num_spectra
    # test file is in the mapping for spectrum "1"
    assert filename in loader.mappings["1"]

    if workflow == GNPSFormat.FBMN:
        assert len(loader.mappings["1"]) == 110
        assert "5425_5426_mod.mzXML" not in loader.mappings["1"]


def test_mapping_reversed(gnps_file_mappings_files):
    loader = GNPSFileMappingLoader(gnps_file_mappings_files[GNPSFormat.SNETSV2])
    assert len(loader.mapping_reversed) == 6
    assert len(loader.mapping_reversed["140221_ME_14_13.mzML"]) == 1028
    assert "1" in loader.mapping_reversed["140221_ME_14_13.mzML"]
