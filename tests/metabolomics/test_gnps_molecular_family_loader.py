import pytest
from nplinker.metabolomics.gnps import GNPSFormat
from nplinker.metabolomics.gnps import GNPSMolecularFamilyLoader


@pytest.mark.parametrize(
    "workflow, num_families, num_spectra",
    [(GNPSFormat.SNETS, 25769, 19), (GNPSFormat.SNETSV2, 6902, 10), (GNPSFormat.FBMN, 1105, 5)],
)
def test_has_molecular_families(workflow, num_families, num_spectra, gnps_mf_files):
    loader = GNPSMolecularFamilyLoader(gnps_mf_files[workflow])
    actual = loader.families
    assert len(actual) == num_families
    # test molecular family with id "1" has correct number of spectra ids
    mf = [mf for mf in actual if mf.family_id == "1"][0]
    assert len(mf.spectra_ids) == num_spectra
