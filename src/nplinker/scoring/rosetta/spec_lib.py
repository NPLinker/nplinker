# Copyright 2021 The NPLinker Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import logging
from sortedcontainers import SortedList
from nplinker.metabolomics.gnps import GNPSSpectrumLoader
from .rosetta_functions import fast_cosine


logger = logging.getLogger(__name__)


class SpecLib:
    def __init__(self, mgf_file):
        self.mgf_file = mgf_file
        self.spectra = []

    def _load_mgf(self):
        self.spectra = GNPSSpectrumLoader(self.mgf_file).spectra()

        self.sort()

    def sort(self):
        # make a sorted list for quick precursor matching
        self.sorted_spectra = [s for s in self.spectra]
        self.sorted_spectra.sort()

    def get_n_spec(self):
        return len(self.spectra)

    def get_ids(self):
        return list(s.id for s in self.spectra)

    def get_n_peaks(self):
        return [len(s.peaks) for s in self.spectra]

    def filter(self):
        # top_k_filter
        n_done = 0
        for spec in self.spectra:
            self._keep_top_k(spec)
            n_done += 1
            if n_done % 100 == 0:
                logger.info(
                    "SpecLib filtered {}/{}, {:.2f}%".format(
                        n_done, len(self.spectra), 100 * (n_done / len(self.spectra))
                    )
                )

    def spectral_match(
        self,
        query,
        scoring_function=fast_cosine,
        ms2_tol=0.2,
        min_match_peaks=1,
        ms1_tol=0.2,
        score_thresh=0.7,
    ):
        candidates = self._candidates(self.sorted_spectra, query.precursor_mz, ms1_tol)
        hits = []
        for c in candidates:
            sc, _ = scoring_function(query, c, ms2_tol, min_match_peaks)
            if sc >= score_thresh:
                hits.append((c.gnps_id, sc))
        return hits

    def _candidates(self, mz_list, query_mz, ms1_tol):
        pmz_list = SortedList([m.precursor_mz for m in mz_list])
        lower = query_mz - ms1_tol
        upper = query_mz + ms1_tol
        start = pmz_list.bisect(lower)
        end = pmz_list.bisect(upper)
        return mz_list[start:end]

    # from molnet repo
    def _keep_top_k(self, spec, k=6, mz_range=50):
        # only keep peaks that are in the top k in += mz_range
        start_pos = 0
        new_mz = []
        new_intensities = []
        for mz, intensity in spec.peaks:
            while spec.peaks[start_pos][0] < mz - mz_range:
                start_pos += 1
            end_pos = start_pos

            n_bigger = 0
            while end_pos < len(spec.peaks) and spec.peaks[end_pos][0] <= mz + mz_range:
                if spec.peaks[end_pos][1] > intensity:
                    n_bigger += 1
                end_pos += 1

            if n_bigger < k:
                new_mz.append(mz)
                new_intensities.append(intensity)

        spec.mz = new_mz
        spec.intensity = new_intensities
