import pytest
from nplinker.metabolomics.gnps import GNPSFormat
from nplinker.metabolomics.gnps import GNPSMolecularFamilyLoader


@pytest.mark.parametrize(
    "workflow, num_families, num_spectra, keep_singleton",
    [
        (GNPSFormat.SNETS, 25769, 19, True),
        (GNPSFormat.SNETSV2, 6902, 10, True),
        (GNPSFormat.FBMN, 1105, 5, True),
        (GNPSFormat.SNETS, 29, 19, False),
        (GNPSFormat.SNETSV2, 72, 10, False),
        (GNPSFormat.FBMN, 60, 5, False),
    ],
)
def test_gnps_molecular_family_loader(
    workflow, num_families, num_spectra, keep_singleton, gnps_mf_files
):
    """Test GNPSMolecularFamilyLoader class."""
    loader = GNPSMolecularFamilyLoader(gnps_mf_files[workflow])
    actual = loader.get_mfs(keep_singleton=keep_singleton)
    assert len(actual) == num_families
    # test molecular family with id "1" has correct number of spectra ids
    mf = [mf for mf in actual if mf.family_id == "1"][0]
    assert len(mf.spectra_ids) == num_spectra
