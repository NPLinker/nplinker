import re

SEARCH_OPT_BGC_NAME,\
SEARCH_OPT_BGC_STRAIN,\
SEARCH_OPT_BGC_BIGSCAPE_CLASS,\
SEARCH_OPT_BGC_PRODUCT_PREDICTION,\
SEARCH_OPT_BGC_DESCRIPTION,\
SEARCH_OPT_GCF_ID,\
SEARCH_OPT_GCF_STRAINS,\
SEARCH_OPT_SPEC_ID,\
SEARCH_OPT_SPEC_ANNOTATIONS,\
SEARCH_OPT_SPEC_FAMILY,\
SEARCH_OPT_SPEC_STRAINS,\
SEARCH_OPT_MOLFAM_FAMILY,\
SEARCH_OPT_MOLFAM_SPECTRA = range(13)

SEARCH_OPTIONS = [
    'BGC-Name',
    'BGC-Strain',
    'BGC-BigscapeClass',
    'BGC-ProductPrediction',
    'BGC-Description',
    'GCF-ID',
    'GCF-Strains',
    'Spectrum-ID',
    'Spectrum-Annotations',
    'Spectrum-Family',
    'Spectrum-Strains',
    'MolFam-Family',
    'MolFam-Spectra',
]

class Searcher(object):

    def __init__(self, nplinker):
        self.npl = nplinker
        self.results = []

    def name_to_mode(self, name):
        return SEARCH_OPTIONS.index(name)

    def textmatch(self, needle, haystack, use_re, exact=False):
        if haystack is None:
            return False

        if use_re:
            return (needle.search(haystack) is not None)

        if exact:
            return (needle == haystack)

        return (haystack.find(needle) != -1)

    def textmatch_list(self, needle, haystacks, use_re):
        for haystack in haystacks:
            if self.textmatch(needle, haystack, use_re):
                return True

        return False

    def search(self, mode, text, use_re=True):
        # TODO all modes
        results = []
        needle = text if not use_re else re.compile(text)

        mode = self.name_to_mode(mode)
        if mode == SEARCH_OPT_BGC_STRAIN:
            results = [bgc for bgc in self.npl.bgcs if self.textmatch(needle, bgc.strain, use_re)]
        elif mode == SEARCH_OPT_BGC_NAME:
            results = [bgc for bgc in self.npl.bgcs if self.textmatch(needle, bgc.name, use_re)]
        elif mode == SEARCH_OPT_BGC_BIGSCAPE_CLASS:
            results = [bgc for bgc in self.npl.bgcs if self.textmatch(needle, bgc.bigscape_class, use_re)]
        elif mode == SEARCH_OPT_BGC_PRODUCT_PREDICTION:
            results = [bgc for bgc in self.npl.bgcs if self.textmatch(needle, ' '.join(bgc.product_prediction), use_re)]
        elif mode == SEARCH_OPT_BGC_DESCRIPTION:
            results = [bgc for bgc in self.npl.bgcs if self.textmatch(needle, bgc.description, use_re)]
        elif mode == SEARCH_OPT_GCF_STRAINS:
            results = [gcf for gcf in self.npl.gcfs if self.textmatch_list(needle, [bgc.strain for bgc in gcf.bgc_list], use_re)]
        elif mode == SEARCH_OPT_GCF_ID:
            results = [gcf for gcf in self.npl.gcfs if self.textmatch(needle, gcf.gcf_id, use_re)]
        elif mode == SEARCH_OPT_SPEC_STRAINS:
            results = [spec for spec in self.npl.spectra if self.textmatch_list(needle, spec.strain_list, use_re)]
        elif mode == SEARCH_OPT_SPEC_FAMILY:
            # if not in regex mode here, probably want to use exact matching
            results = [spec for spec in self.npl.spectra if self.textmatch(needle, spec.family, use_re, not use_re)]
        elif mode == SEARCH_OPT_SPEC_ID:
            # if not in regex mode here, probably want to use exact matching
            results = [spec for spec in self.npl.spectra if self.textmatch(needle, str(spec.id), use_re, not use_re)]

        print('found {} results'.format(len(results)))
        self.results = results
        return results
