import pytest
from nplinker.metabolomics.gnps import GNPSFormat
from nplinker.metabolomics.gnps import GNPSSpectrumLoader


@pytest.mark.parametrize("workflow, num_spectra",
                         [[GNPSFormat.FBMN, 1492], [GNPSFormat.SNETS, 25935],
                          [GNPSFormat.SNETSV2, 7383]])
def test_gnps_spectrum_loader(workflow, num_spectra, gnps_spectra_files):
    loader = GNPSSpectrumLoader(gnps_spectra_files[workflow])
    assert len(loader.spectra) == num_spectra
