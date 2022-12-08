from os import PathLike
from nplinker.annotations import GNPS_DATA_COLUMNS
from nplinker.annotations import GNPS_KEY
from nplinker.annotations import create_gnps_annotation
from nplinker.logconfig import LogConfig
from nplinker.metabolomics.abc import SpectrumLoaderBase
from nplinker.metabolomics.spectrum import Spectrum
from nplinker.parsers.mgf import LoadMGF


logger = LogConfig.getLogger(__file__)

class GNPSSpectrumLoader(SpectrumLoaderBase):

    def __init__(self, filename: str | PathLike):
        """Load the mass spectra from the MGF file pointed to by `filename`.

        Args:
            filename(str | PathLike): str or PathLike object pointing to the spectra files to load.
        """
        ms1, ms2, metadata = LoadMGF(name_field='scans').load_spectra([str(filename)])
        logger.info('%d molecules parsed from MGF file', len(ms1))
        self._spectra = _mols_to_spectra(ms2, metadata)
    
    def spectra(self) -> list[Spectrum]:
        """Get the spectra loaded from the file.

        Returns:
            list[Spectrum]: Spectra loaded from the file.
        """
        return self._spectra
    

def _mols_to_spectra(ms2, metadata: dict[str, dict[str, str]]) -> list[Spectrum]:
    """Function to convert ms2 object and metadata to `Spectrum` objects.

    Args:
        ms2(_type_): Unknown.
        metadata(dict[str, dict[str, str]]): Dictionary holding the metadata which was loaded from the original file.

    Returns:
        list[Spectrum]: List of mass spectra obtained from ms2 and metadata.
    """
    ms2_dict = {}
    for m in ms2:
        if not m[3] in ms2_dict:
            ms2_dict[m[3]] = []
        ms2_dict[m[3]].append((m[0], m[2]))

    spectra = []
    for i, m in enumerate(ms2_dict.keys()):
        new_spectrum = Spectrum(i, ms2_dict[m], int(m.name),
                                metadata[m.name]['precursormass'],
                                metadata[m.name]['parentmass'])
        new_spectrum.metadata = metadata[m.name]
        # add GNPS ID if in metadata under SPECTRUMID (this doesn't seem to be in regular MGF files
        # from GNPS, but *is* in the rosetta mibig MGF)
        # note: LoadMGF seems to lowercase (some) metadata keys?
        if 'spectrumid' in new_spectrum.metadata:
            # add an annotation for consistency, if not already there
            if GNPS_KEY not in new_spectrum.annotations:
                gnps_anno = {k: None for k in GNPS_DATA_COLUMNS}
                gnps_anno['SpectrumID'] = new_spectrum.metadata['spectrumid']
                create_gnps_annotation(new_spectrum, gnps_anno)
        spectra.append(new_spectrum)

    return spectra
