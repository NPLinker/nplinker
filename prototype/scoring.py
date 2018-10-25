
def compute_all_scores(spectra_list,gcf_list,strain_list,scoring_function,do_random = True):
    m_scores = {}
    best = 0
    best_random = 0
    for i,spectrum in enumerate(spectra_list):
        m_scores[spectrum] = {}
        if i % 100 == 0:
            print "Done {} of {}".format(i,len(spectra_list))
        for gcf in gcf_list:
            s = scoring_function(spectrum,gcf,strain_list)
            if do_random:
                s_random = scoring_function(spectrum.random_spectrum,gcf.random_gcf,strain_list)
            else:
                s_random = None
            m_scores[spectrum][gcf] = (s,s_random,None)
            if s > best:
                best = s
                print "Best: ",best
            if s_random > best_random:
                best_random = s_random
                print "Best random: ",best_random
    return m_scores

def metcalf_scoring(spectral_like,gcf_like,strains,both = 10,met_not_gcf = -10,gcf_not_met = 0,neither = 1):
    cum_score = 0
    for strain in strains:
        in_spec = spectral_like.has_strain(strain)
        in_gcf = gcf_like.has_strain(strain)
        if in_spec and in_gcf:
            cum_score += both
        if in_spec and not in_gcf:
            cum_score += met_not_gcf
        if in_gcf and not in_spec:
            cum_score += gcf_not_met
        if not in_gcf and not in_spec:
            cum_score += 1
    return cum_score

def name_scoring(spectral_like,gcf_like,strains):
	if len(spectral_like.annotations) == 0:
		return None
	spectral_annotations = set()
	for annotation in spectral_like.annotations:
		spectral_annotations.add(annotation)
	mibig_bgcs = []
	for bgc in gcf_like.bgc_list:
		print type(bgc)
