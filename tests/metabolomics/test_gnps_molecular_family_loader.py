import os
import pytest
from nplinker.metabolomics.gnps.gnps_molecular_family_loader import \
    GNPSMolecularFamilyLoader
from .. import DATA_DIR


@pytest.mark.parametrize(
    "filename",
    [os.path.join(DATA_DIR, "edges.pairsinfo"), DATA_DIR / "edges.pairsinfo"])
def test_has_molecular_families(filename):
    sut = GNPSMolecularFamilyLoader(filename)
    actual = sut.families()
    assert len(actual) == 25769
    mf_ids = [mf.family_id for mf in actual[:30]]
    assert mf_ids == [
        '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13',
        '14', '15', '16', '17', '18', '20', '21', '22', '23', '24', '26', '28',
        '30', '31', '32', '33', '-1'
    ]
    num_spec_ids = [len(mf.spectra_ids) for mf in actual[:30]]
    assert num_spec_ids == [
        19, 48, 3, 3, 11, 4, 9, 3, 15, 3, 5, 2, 3, 3, 5, 3, 14, 4, 2, 2, 12, 2,
        3, 5, 2, 4, 2, 2, 2, 1
    ]
    assert actual[0].spectra_ids == set(
        ('13170', '13662', '15316', '15364', '16341', '17201', '17270',
         '18120', '18172', '18748', '18831', '19005', '19673', '19719',
         '20320', '20738', '21163', '21593', '23566'))
