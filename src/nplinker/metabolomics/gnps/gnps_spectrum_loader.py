from os import PathLike
from nplinker.logconfig import LogConfig
from nplinker.metabolomics.abc import SpectrumLoaderBase
from nplinker.metabolomics.spectrum import Spectrum
from nplinker.parsers.mgf import LoadMGF


logger = LogConfig.getLogger(__file__)

class GNPSSpectrumLoader(SpectrumLoaderBase):

    def __init__(self, file: str | PathLike):
        """Class to load mass spectra from the given GNPS MGF file.

        Args:
            file(str | PathLike): path to the MGF file to load.
        """
        ms1, ms2, metadata = LoadMGF(name_field='scans').load_spectra([str(file)])
        logger.info('%d molecules parsed from MGF file', len(ms1))
        self._spectra = _mols_to_spectra(ms2, metadata)

    def spectra(self) -> list[Spectrum]:
        """Get the spectra loaded from the file.

        Returns:
            list[Spectrum]: the loaded spectra as a list of `Spectrum` objects.
        """
        return self._spectra


def _mols_to_spectra(ms2: list, metadata: dict[str, dict[str, str]]) -> list[Spectrum]:
    """Function to convert ms2 object and metadata to `Spectrum` objects.

    Args:
        ms2(list): list of ms2 intensity.
        metadata(dict[str, dict[str, str]]): Dictionary holding the metadata which was loaded from the original file.

    Returns:
        list[Spectrum]: List of mass spectra obtained from ms2 and metadata.
    """
    ms2_dict = {}
    # an example of m:
    # (118.487999, 0.0, 18.753, <nplinker.parsers.mg...105f2c970>, 'spectra.mgf', 0.0)
    for m in ms2:
        if not m[3] in ms2_dict: # m[3] is `nplinker.parsers.mgf.MS1` object
            ms2_dict[m[3]] = []
        ms2_dict[m[3]].append((m[0], m[2]))

    spectra = []
    for i, m in enumerate(ms2_dict.keys()): # m is `nplinker.parsers.mgf.MS1` object
        new_spectrum = Spectrum(i, ms2_dict[m], m.name,
                                metadata[m.name]['precursormass'],
                                metadata[m.name]['parentmass'])
        new_spectrum.metadata = metadata[m.name]
        spectra.append(new_spectrum)

    return spectra
