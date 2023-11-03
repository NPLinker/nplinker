from nplinker.strain_collection import StrainCollection
from nplinker.strains import Strain
from nplinker.utils import sqrt_normalise


GNPS_KEY = "gnps"

JCAMP = (
    "##TITLE={}\\n"
    + "##JCAMP-DX=nplinker vTODO\\n"
    + "##DATA TYPE=Spectrum\\n"
    + "##DATA CLASS=PEAKTABLE\\n"
    + "##ORIGIN=TODO_DATASET_ID\\n"
    + "##OWNER=nobody\\n"
    + "##XUNITS=M/Z\\n"
    + "##YUNITS=RELATIVE ABUNDANCE\\n"
    + "##NPOINTS={}\\n"
    + "##PEAK TABLE=(XY..XY)\\n"
    + "{}\\n"
    + "##END=\\n"
)


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
        self.growth_media = {}
        # TODO CG: self.family_id should be removed, used in deprecated make_families method
        self.family_id = "-1"
        self.family = None
        # a dict indexed by filename, or "gnps"
        self.annotations = {}
        self._losses = None
        self._jcamp = None

    def add_strain(self, strain, growth_medium, peak_intensity):
        # adds the strain to the StrainCollection if not already there
        self.strains.add(strain)

        if strain not in self.growth_media:
            self.growth_media[strain] = {}

        if growth_medium is None:
            self.growth_media[strain].update(
                {f"unknown_medium_{len(self.growth_media[strain])}": peak_intensity}
            )
            return

        if strain in self.growth_media and growth_medium in self.growth_media[strain]:
            raise Exception("Growth medium clash: {} / {} {}".format(self, strain, growth_medium))

        self.growth_media[strain].update({growth_medium: peak_intensity})

    @property
    def is_library(self):
        return GNPS_KEY in self.annotations

    def set_annotations(self, key, data):
        self.annotations[key] = data

    @property
    def gnps_annotations(self):
        if GNPS_KEY not in self.annotations:
            return None

        return self.annotations[GNPS_KEY][0]

    def has_annotations(self):
        return len(self.annotations) > 0

    def get_metadata_value(self, key):
        val = self.metadata.get(key, None)
        return val

    def has_strain(self, strain: Strain):
        return strain in self.strains

    def get_growth_medium(self, strain):
        if strain not in self.strains:
            return None

        gms = self.growth_media[strain]
        return list(gms.keys())[0]

    def to_jcamp_str(self, force_refresh=False):
        if self._jcamp is not None and not force_refresh:
            return self._jcamp

        peakdata = "\\n".join("{}, {}".format(*p) for p in self.peaks)
        self._jcamp = JCAMP.format(str(self), self.n_peaks, peakdata)
        return self._jcamp

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

    @property
    def losses(self):
        """All mass shifts in the spectrum, and the indices of the peaks."""
        if self._losses is None:
            # populate loss table
            losses = []
            for i in range(len(self.peaks)):
                loss = self.precursor_mz - self.peaks[i][0]
                losses.append((loss, self.id, i))

            # THIS SEEMED TO ME LIKE IT WOULD TAKE THE WRONG DIFFERENCES AS LOSSES:
            # TODO: please check!
            #                for j in range(i):
            #                    loss = self.peaks[i][0] - self.peaks[j][0]
            #                    losses.append((loss, i, j))

            # Sort by loss
            losses.sort(key=lambda x: x[0])
            self._losses = losses
        return self._losses

    def has_loss(self, mass, tol):
        """Check if the scan has the specified loss (within tolerance)."""
        matched_losses = []

        idx = 0
        # Check losses in range [0, mass]
        while idx < len(self.losses) and self.losses[idx][0] <= mass:
            if mass - self.losses[idx][0] < tol:
                matched_losses.append(self.losses[idx])
            idx += 1

        # Add all losses in range [mass, mass+tol(
        while idx < len(self.losses) and self.losses[idx][0] < mass + tol:
            matched_losses.append(self.losses[idx])
            idx += 1

        return matched_losses
