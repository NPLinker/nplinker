# some code for spectral library things

from .scoring_functions import fast_cosine, fast_cosine_shift

class SpecLib(object):
    def __init__(self,mgf_file):
        self.mgf_file = mgf_file
        self.spectra = None

    def _load_mgf(self,id_field='SPECTRUMID'):
        from mnet_utilities import load_mgf
        self.spectra = load_mgf(self.mgf_file,id_field = id_field)
        for k,v in self.spectra.items():
            v.spectrum_id = k

    def get_n_spec(self):
        return len(self.spectra)

    def get_keys(self):
        return list(self.spectra.keys())

    def get_n_peaks(self):
        return [s.n_peaks for s in self.spectra.values()]

    def filter(self):
        # top_k_filter
        n_done = 0
        for s_id,spec in self.spectra.items():
            spec.keep_top_k()
            n_done += 1
            if n_done % 100 == 0:
                print("Filtered {}".format(n_done))

    
    def spectral_match(self,query,
            scoring_function = fast_cosine,
            ms2_tol = 0.2,
            min_match_peaks = 1,
            ms1_tol = 0.2,
            score_thresh = 0.7):
        # make a sorted list for quick precursor matching
        spec_list = [s for s in self.spectra.values()]
        spec_list.sort()
        candidates = self._candidates(spec_list,query.precursor_mz,ms1_tol)
        hits = []
        for c in candidates:
            sc,_ = scoring_function(query,c,ms2_tol,min_match_peaks)
            if sc >= score_thresh:
                hits.append((c.spectrum_id,sc))
        return hits

        
    def _candidates(self,mz_list,query_mz,ms1_tol):
        from sortedcontainers import SortedList
        pmz_list = SortedList([m.precursor_mz for m in mz_list])
        lower = query_mz - ms1_tol
        upper = query_mz + ms1_tol
        start = pmz_list.bisect(lower)
        end = pmz_list.bisect(upper)
        return mz_list[start:end]
        

