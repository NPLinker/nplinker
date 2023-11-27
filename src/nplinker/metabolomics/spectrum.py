from __future__ import annotations
from typing import TYPE_CHECKING
from nplinker.strain import Strain
from nplinker.strain_collection import StrainCollection
from nplinker.utils import sqrt_normalise


if TYPE_CHECKING:
    from .molecular_family import MolecularFamily

GNPS_KEY = "gnps"


class Spectrum:
    def __init__(self, id, peaks, spectrum_id: str, precursor_mz, parent_mz=None, rt=None):
        self.id = id
        self.peaks = sorted(peaks, key=lambda x: x[0])  # ensure sorted by mz
        self.normalised_peaks = sqrt_normalise(self.peaks)  # useful later
        self.n_peaks = len(self.peaks)
        self.max_ms2_intensity = max(intensity for mz, intensity in self.peaks)
        self.total_ms2_intensity = sum(intensity for mz, intensity in self.peaks)
        self.spectrum_id = spectrum_id  # MS1.name
        self.rt = rt
        # TODO CG: should include precursor mass and charge to calculate precursor_mz
        # parent_mz can be calculate from precursor_mass and charge mass
        self.precursor_mz = precursor_mz
        self.parent_mz = parent_mz
        self.gnps_id = None  # CCMSLIB...
        # TODO should add intensity here too
        self.metadata = {}
        self.edges = []
        self.strains = StrainCollection()
        # this is a dict indexed by Strain objects (the strains found in this Spectrum), with
        # the values being dicts of the form {growth_medium: peak intensity} for the parent strain
        self.family: MolecularFamily | None = None
        # a dict indexed by filename, or "gnps"
        self.annotations = {}

    @property
    def gnps_annotations(self):
        if GNPS_KEY not in self.annotations:
            return None

        return self.annotations[GNPS_KEY][0]

    def has_strain(self, strain: Strain):
        return strain in self.strains

    def __str__(self):
        return "Spectrum(id={}, spectrum_id={}, strains={})".format(
            self.id, self.spectrum_id, len(self.strains)
        )

    def __repr__(self):
        return str(self)

    def __eq__(self, other) -> bool:
        if isinstance(other, Spectrum):
            return (
                self.id == other.id
                and self.spectrum_id == other.spectrum_id
                and self.precursor_mz == other.precursor_mz
                and self.parent_mz == other.parent_mz
            )
        return NotImplemented

    def __hash__(self) -> int:
        return hash((self.id, self.spectrum_id, self.precursor_mz, self.parent_mz))

    def __cmp__(self, other):
        if self.parent_mz >= other.parent_mz:
            return 1
        else:
            return -1

    def __lt__(self, other):
        if self.parent_mz <= other.parent_mz:
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
        self.n_peaks = len(self.peaks)
        if self.n_peaks > 0:
            self.normalised_peaks = sqrt_normalise(self.peaks)
            self.max_ms2_intensity = max(intensity for mz, intensity in self.peaks)
            self.total_ms2_intensity = sum(intensity for mz, intensity in self.peaks)
        else:
            self.normalised_peaks = []
            self.max_ms2_intensity = 0.0
            self.total_ms2_intensity = 0.0
