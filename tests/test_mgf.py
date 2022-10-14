import os
from nplinker.parsers.mgf import LoadMGF


def test_load_mgf():
    loader = LoadMGF(name_field='scans')
    mgf_filepath = os.path.join(os.getcwd(),"tests", "data", "METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra-main.mgf")
    ms1, ms2, metadata = loader.load_spectra([mgf_filepath])

    assert ms1 is not None
    assert ms2 is not None
    assert metadata is not None