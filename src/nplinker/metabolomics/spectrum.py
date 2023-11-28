from __future__ import annotations
from typing import TYPE_CHECKING
from nplinker.strain import Strain
from nplinker.strain_collection import StrainCollection
from nplinker.utils import sqrt_normalise


if TYPE_CHECKING:
    from .molecular_family import MolecularFamily

GNPS_KEY = "gnps"


class Spectrum:
    def __init__(
        self,
        spectrum_id: str,
        peaks: list[tuple[float, float]],
        precursor_mz: float,
        rt: float | None = None,
        metadata: dict | None = None,
        annotations: dict | None = None,
    ) -> None:
        """Class to model MS/MS Spectrum.

        Args:
            spectrum_id (str): the spectrum ID.
            peaks (list[tuple[float, float]]): the list of ordered peaks, as tuples of
                (m/z, intensity). Make sure the peaks are sorted by m/z.
            precursor_mz (float): the precursor m/z.
            rt (float, optional): the retention time in seconds.
            metadata (dict, optional): the metadata of the spectrum, i.e. the header infomation
                in the MGF file.
            annotations (dict, optional): the annotations of the spectrum, e.g. annotations from
                GNPS.

        Attributes:
            spectrum_id (str): the spectrum ID.
            peaks (list[tuple[float, float]]): the list of ordered peaks, as tuples of
                (m/z, intensity), ordered by m/z.
            precursor_mz (float): the precursor m/z.
            rt (float): the retention time in seconds.
            metadata (dict): the metadata of the spectrum, i.e. the header infomation in the MGF
                file.
            annotations (dict): the annotations of the spectrum, e.g. annotations from GNPS.
            normalised_peaks (list[tuple[float, float]]): the list of normalised peaks, ordered by
                m/z.
            gnps_id (str): the GNPS ID of the spectrum.
            strains (StrainCollection): the strains that this spectrum belongs to.
            family (MolecularFamily): the molecular family that this spectrum belongs to.
        """
        self.spectrum_id = spectrum_id
        self.peaks = peaks
        self.precursor_mz = precursor_mz
        self.rt = rt
        self.metadata = metadata or {}
        self.annotations = annotations or {}

        self.normalised_peaks = sqrt_normalise(self.peaks)

        self.gnps_id = None
        self.strains = StrainCollection()
        self.family: MolecularFamily | None = None

    @property
    def gnps_annotations(self):
        if GNPS_KEY not in self.annotations:
            return None

        return self.annotations[GNPS_KEY][0]

    def has_strain(self, strain: Strain):
        return strain in self.strains

    def __str__(self):
        return "Spectrum(spectrum_id={}, strains={})".format(self.spectrum_id, len(self.strains))

    def __repr__(self):
        return str(self)

    def __eq__(self, other) -> bool:
        if isinstance(other, Spectrum):
            return self.spectrum_id == other.spectrum_id and self.precursor_mz == other.precursor_mz
        return NotImplemented

    def __hash__(self) -> int:
        return hash((self.spectrum_id, self.precursor_mz))

    def __cmp__(self, other):
        if self.precursor_mz >= other.precursor_mz:
            return 1
        else:
            return -1

    def __lt__(self, other):
        if self.precursor_mz <= other.precursor_mz:
            return 1
        else:
            return 0

    # from molnet repo
    def keep_top_k(self, k=6, mz_range=50):
        # only keep peaks that are in the top k in += mz_range
        start_pos = 0
        new_peaks = []
        for mz, intensity in self.peaks:
            while self.peaks[start_pos][0] < mz - mz_range:
                start_pos += 1
            end_pos = start_pos

            n_bigger = 0
            while end_pos < len(self.peaks) and self.peaks[end_pos][0] <= mz + mz_range:
                if self.peaks[end_pos][1] > intensity:
                    n_bigger += 1
                end_pos += 1

            if n_bigger < k:
                new_peaks.append((mz, intensity))

        self.peaks = new_peaks
        if len(self.peaks) > 0:
            self.normalised_peaks = sqrt_normalise(self.peaks)
        else:
            self.normalised_peaks = []
