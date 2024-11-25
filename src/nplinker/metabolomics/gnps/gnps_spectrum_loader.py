from __future__ import annotations
import logging
from os import PathLike
from pyteomics import mgf
from nplinker.metabolomics import Spectrum
from nplinker.metabolomics.abc import SpectrumLoaderBase


logger = logging.getLogger(__name__)


class GNPSSpectrumLoader(SpectrumLoaderBase):
    """Load mass spectra from the given GNPS MGF file.

    ??? info "Concept"
        [GNPS data][gnps-data]

    The file mappings file is from GNPS output archive, as described below
    for each GNPS workflow type:

    1. METABOLOMICS-SNETS
        - METABOLOMICS-SNETS*.mgf
    2. METABOLOMICS-SNETS-V2
        - METABOLOMICS-SNETS-V2*.mgf
    3. FEATURE-BASED-MOLECULAR-NETWORKING
        - spectra/*.mgf
    """

    def __init__(self, file: str | PathLike) -> None:
        """Initialize the GNPSSpectrumLoader.

        Args:
            file: path to the MGF file.

        Raises:
            ValueError: Raises ValueError if the file is not valid.

        Examples:
            >>> loader = GNPSSpectrumLoader("gnps_spectra.mgf")
            >>> print(loader.spectra[0])
        """
        self._file = str(file)
        self._spectra: list[Spectrum] = []

        self._validate()
        self._load()

    @property
    def spectra(self) -> list[Spectrum]:
        """Get the list of Spectrum objects.

        Returns:
            list[Spectrum]: the loaded spectra as a list of `Spectrum` objects.
        """
        return self._spectra

    def _validate(self) -> None:
        """Validate GNPS MGF file.

        Raises:
            ValueError: Raises ValueError if the file is not valid.
        """
        # check the local scope of a single MS/MS query (spectrum) has the
        # required parameters. Note that this is not the header of the MGF
        # file, but the local scope of each spectrum.
        required_params = ["scans", "pepmass", "charge"]
        for spec in mgf.MGF(self._file):
            for param in required_params:
                if param not in spec["params"]:
                    raise ValueError(
                        f"Invalid MGF file '{self._file}'. "
                        f"Expected parameter '{param}' not found, "
                        f"but got '{spec['params']}'."
                    )

    def _load(self) -> None:
        """Load the MGF file into Spectrum objects."""
        for spec in mgf.MGF(self._file):
            # Skip if m/z array is empty, as this is an invalid spectrum.
            # The invalid spectrum does not exist in other GNPS files, e.g.
            # file mappings file and molecular families file. So we can safely
            # skip it.
            if len(spec["m/z array"]) == 0:
                continue

            # Load the spectrum
            spectrum_id: str = spec["params"]["scans"]
            # The pepmass in an mgf file is actually the m/z and not the peptide mass 
            # See: https://www.matrixscience.com/help/obsolete_data_file_formats.html
            precursor_mz: float = spec["params"]["pepmass"][0]
            precursor_charge: int = spec["params"]["charge"][0]
            rt = spec["params"].get("rtinseconds", 0)

            spectrum = Spectrum(
                id=spectrum_id,
                mz=list(spec["m/z array"]),
                intensity=list(spec["intensity array"]),
                precursor_mz=precursor_mz,
                precursor_charge=precursor_charge,
                rt=rt,
                metadata=spec["params"],
            )
            self._spectra.append(spectrum)
