from nplinker.annotations import GNPS_DATA_COLUMNS
from nplinker.annotations import GNPS_KEY
from nplinker.annotations import create_gnps_annotation
from nplinker.logconfig import LogConfig
from nplinker.metabolomics.ISpectrumLoader import SpectrumLoaderBase
from nplinker.metabolomics.spectrum import Spectrum
from nplinker.parsers.mgf import LoadMGF


logger = LogConfig.getLogger(__file__)

class GNPSSpectrumLoader(SpectrumLoaderBase):

    def __init__(self, filename: str):
        ms1, ms2, metadata = LoadMGF(name_field='scans').load_spectra([filename])
        logger.info('%d molecules parsed from MGF file', len(ms1))
        self._spectra = mols_to_spectra(ms2, metadata)
    
    def spectra(self) -> list[Spectrum]:
        return self._spectra
    

def mols_to_spectra(ms2, metadata: dict[str, str]) -> list[Spectrum]:
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
