from sortedcontainers import SortedList

from .rosetta_functions import fast_cosine, fast_cosine_shift
from ...parsers.mgf import LoadMGF
from ...metabolomics import mols_to_spectra

from ...logconfig import LogConfig
logger = LogConfig.getLogger(__file__)

class SpecLib(object):
    def __init__(self, mgf_file):
        self.mgf_file = mgf_file
        self.spectra = []

    def _load_mgf(self):
        ms1, ms2, metadata = LoadMGF(name_field='scans').load_spectra([self.mgf_file])
        self.spectra = mols_to_spectra(ms2, metadata)
        self.sort()

    def sort(self):
        # make a sorted list for quick precursor matching
        self.sorted_spectra = [s for s in self.spectra]
        self.sorted_spectra.sort()

    def get_n_spec(self):
        return len(self.spectra)

    def get_ids(self):
        return list(s.spectrum_id for s in self.spectra)

    def get_n_peaks(self):
        return [s.n_peaks for s in self.spectra]

    def filter(self):
        # top_k_filter
        n_done = 0
        for spec in self.spectra:
            spec.keep_top_k()
            n_done += 1
            if n_done % 100 == 0:
                logger.info('SpecLib filtered {}/{}, {:.2f}%'.format(n_done, len(self.spectra), 100 * (n_done / len(self.spectra))))

    def spectral_match(self, query, scoring_function=fast_cosine, ms2_tol=0.2, min_match_peaks=1, ms1_tol=0.2, score_thresh=0.7):
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
        

