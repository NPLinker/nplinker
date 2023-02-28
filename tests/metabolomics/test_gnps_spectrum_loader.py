import pytest
from nplinker.metabolomics.gnps.gnps_spectrum_loader import GNPSSpectrumLoader

from .. import DATA_DIR

@pytest.mark.parametrize("file, expected", [
    [DATA_DIR / "ProteoSAFe-METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra/METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra-main.mgf", 435],
    [DATA_DIR / "ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-92036537-download_cytoscape_data/spectra/specs_ms.mgf", 1492]
])
def test_loads_spectra(file, expected):
    actual = GNPSSpectrumLoader(file).spectra()
    assert len(actual) == expected
