from nplinker.parsers.mgf import LoadMGF
from .. import GNPS_DATA_DIR


def test_load_mgf():
    loader = LoadMGF(name_field="scans")
    ms1, ms2, metadata = loader.load_spectra([GNPS_DATA_DIR / "spectra.mgf"])

    assert ms1 is not None
    assert ms2 is not None
    assert metadata is not None
