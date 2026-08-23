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
def test_file_mapping_loader_gnps1(workflow, num_spectra, filename, gnps_file_mappings_files):
    loader = GNPSFileMappingLoader(gnps_file_mappings_files[workflow])
    assert len(loader.mappings) == num_spectra
    # test file is in the mapping for spectrum "1"
    assert filename in loader.mappings["1"]

    if workflow == GNPSFormat.FBMN:
        assert len(loader.mappings["1"]) == 110
        assert "5425_5426_mod.mzXML" not in loader.mappings["1"]


@pytest.mark.parametrize(
    "workflow, num_spectra, filename",
    [
        [GNPSFormat.GNPS2CN, 1051, "blk_g10_dora.mzML"],
        [GNPSFormat.GNPS2FBMN, 371, "blk_g10_dora.mzML"],
    ],
)
def test_file_mapping_loader_gnps2(workflow, num_spectra, filename, gnps2_file_mappings_files):
    loader = GNPSFileMappingLoader(gnps2_file_mappings_files[workflow])
    assert len(loader.mappings) == num_spectra
    # test file is in the mapping for spectrum "2"
    assert filename in loader.mappings["2"]


def test_mapping_reversed(gnps_file_mappings_files):
    loader = GNPSFileMappingLoader(gnps_file_mappings_files[GNPSFormat.SNETSV2])
    assert len(loader.mapping_reversed) == 6
    assert len(loader.mapping_reversed["140221_ME_14_13.mzML"]) == 1028
    assert "1" in loader.mapping_reversed["140221_ME_14_13.mzML"]


def test_load_csv_trims_peak_area_suffix(tmp_path):
    """Regression: the ' Peak area' column suffix must be removed as a suffix, not
    as a character set.

    `str.strip(" Peak area")` removes any leading/trailing character in
    {space, P, e, a, k, r}, corrupting file names such as
    'extract1.mzML' -> 'xtract1.mzML'. Digit-leading names (as in the real
    fixtures) happen to survive, which is why the bug went unnoticed.
    """
    csv_file = tmp_path / "quantification_table.csv"
    csv_file.write_text(
        "row ID,extract1.mzML Peak area,5434_mod.mzXML Peak area\n1,100.0,0.0\n2,0.0,55.0\n"
    )
    loader = GNPSFileMappingLoader(csv_file)
    # leading 'e' must survive (str.strip would eat it), and the value==0 column excluded
    assert loader.mappings["1"] == ["extract1.mzML"]
    # digit-leading name unchanged (no regression)
    assert loader.mappings["2"] == ["5434_mod.mzXML"]
