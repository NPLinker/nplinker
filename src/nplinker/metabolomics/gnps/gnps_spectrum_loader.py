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
            # calculate precursor m/z from precursor mass and charge
            precursor_mass = spec["params"]["pepmass"][0]
            precursor_charge = self._get_precursor_charge(spec["params"]["charge"])
            precursor_mz: float = precursor_mass / abs(precursor_charge)
            rt = spec["params"].get("rtinseconds", 0)

            spectrum = Spectrum(
                id=spectrum_id,
                mz=list(spec["m/z array"]),
                intensity=list(spec["intensity array"]),
                precursor_mz=precursor_mz,
                rt=rt,
                metadata=spec["params"],
            )
            self._spectra.append(spectrum)

    def _get_precursor_charge(self, charges: list[int]) -> int:
        """Get the precursor charge from the charge list.

        Args:
            charges: list of charge values.

        Returns:
            the precursor charge.
        """
        charge = charges[0]
        if charge == 0:
            logger.warning(
                f"Invalid precursor charge value 0. "
                f"Assuming charge is 1 for spectrum '{self._file}'."
            )
            charge = 1
        return charge
