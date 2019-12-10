import numpy as np
import scipy.stats

# weights for metcalf_scoring, in order:
# - neither
# - GCF only
# - spectrum only
# - both
METCALF_WEIGHTS = [1, 0, -10, 10]


METCALF_EXPECTED_CACHE = None
METCALF_VARIANCE_CACHE = None


def metcalf_std_count(*arg):
    raw = metcalf_count(*arg)
    exp, std = metcalf_expected_count(*arg)
    return (raw - exp) / np.sqrt(std)


def metcalf_count(n, m, overlap, N, in_both, spec_not_gen, gen_not_spec, neither):
    # compute metcalf score for n strains in spectrum
    # m strains in gcf
    # where overlap of them are the same
    score = overlap * in_both
    score += spec_not_gen*(n-overlap)
    score += gen_not_spec*(m-overlap)
    score += neither * (N - (n + m - overlap))
    return score


def metcalf_expected_count(met, gcf, overlap, N, both, met_not_gcf, gcf_not_met, neither):
    # Compute the expected value of the Metcalf score for the given MF/GCF sizes
    global METCALF_EXPECTED_CACHE
    global METCALF_VARIANCE_CACHE

    cache_size = (N + 1, N + 1)
    if METCALF_EXPECTED_CACHE is None or METCALF_EXPECTED_CACHE.shape != cache_size:
        METCALF_EXPECTED_CACHE = np.zeros(cache_size)
        METCALF_VARIANCE_CACHE = np.zeros(cache_size)

    if METCALF_EXPECTED_CACHE[met, gcf] != 0.0 or METCALF_VARIANCE_CACHE[met, gcf] != 0.0:
        expected = METCALF_EXPECTED_CACHE[met, gcf]
        variance = METCALF_VARIANCE_CACHE[met, gcf]
        return expected, variance

    min_overlap = max(0, met + gcf - N)
    max_overlap = min(met, gcf)
    this_e = 0
    this_e2 = 0
    for o in range(int(min_overlap), int(max_overlap)+1):
        o_prob = scipy.stats.hypergeom.pmf(o, N, met, gcf)
        score = metcalf_count(met, gcf, o, N, both, met_not_gcf, gcf_not_met, neither)
        this_e += o_prob*score
        this_e2 += o_prob*(score**2)
                
    expected = this_e
    variance = this_e2 - this_e**2
    # degenerate example...
    if variance == 0:
        variance = 1

    METCALF_EXPECTED_CACHE[met, gcf] = expected
    METCALF_VARIANCE_CACHE[met, gcf] = variance

    return expected, variance


def metcalf_scoring(spectral_like, gcf_like, strains, both=10, met_not_gcf=-10, gcf_not_met=0, neither=1, standardised=False):
    """
    Calculate the Metcalf score for the provided objects
    strains should be a list of strain names
    """
    N = len(strains)

    mf_strains = list(set(spectral_like.strains).intersection(strains))
    gcf_strains = list(set(gcf_like.strains).intersection(strains))

    n = len(mf_strains)
    m = len(gcf_strains)
    o = len(set(mf_strains).intersection(gcf_strains))
    actual_score = metcalf_count(n, m, o, N, both, met_not_gcf, gcf_not_met, neither)

    if standardised:
        expected, variance = metcalf_expected_count(n, m, o, N, both, met_not_gcf, gcf_not_met, neither)
        actual_score = (actual_score - expected) / np.sqrt(variance)

    return actual_score


def hg_scoring(spectral_like, gcf_like, strains):
    spectral_count = 0
    gcf_count = 0
    overlap_count = 0

    for strain in strains:
        in_spec = spectral_like.has_strain(strain)
        in_gcf = gcf_like.has_strain(strain)
        if in_spec:
            spectral_count += 1
        if in_gcf:
            gcf_count += 1
        if in_gcf and in_spec:
            overlap_count += 1

    pos_in_sample = overlap_count
    N = spectral_count
    n = gcf_count
    M = len(strains)

    # sf(k, M, n, N, loc=0)
    # sf(num_shared_strains, num_strains, strains_in_gcf, strains_in_spec, 1)
    r = scipy.stats.hypergeom.sf(pos_in_sample, M, n, N, 1)
    return r, None

# TODO needs updating due to annotation changes
def name_scoring(spectral_like, gcf_like, mibig_map):
    score = 0
    metadata = None
    if spectral_like.empty_default_annotations():
        print("No annotations")
        return None, None

    mibig_bgcs = gcf_like.mibig_bgcs
    if len(mibig_bgcs) == 0:
        print("no mibig")
        return None, None

    for annotation in spectral_like.get_annotations():
        for mibig in mibig_bgcs:
            short_mibig = mibig.name.split('.')[0]
            if short_mibig in mibig_map:
                m = match(annotation, mibig_map[short_mibig])
                if m:
                    metadata = m
                    score += 100
    return (score, metadata)

def match(spectral_annotation, mibig_name):
    metadata = None
    name, source = spectral_annotation
    for m_name in mibig_name:
        if name.lower() == m_name.split()[0].lower():
            print(name, m_name)
            metadata = (name, m_name)
            return metadata
    return False

# TODO needs updating due to annotation changes
def knownclusterblast_scoring(spectral_like, gcf_like, mibig_map):
    score = 0
    metadata = None
    if spectral_like.empty_default_annotations():
        print("No annotations")
        return None, None
    kcb = []
    for bgc in gcf_like.bgc_list:
        these = bgc.known_cluster_blast
        # if hasattr(bgc,'metadata'):
        #     these = bgc.metadata.get('knownclusterblast',None)
        if these:
            for mibig, score in these:
                kcb.append((mibig, score))
    if len(kcb) == 0:
        return None, None
    total_score = 0
    for annotation in spectral_like.get_annotations():
        for mibig, score in kcb:
            short_mibig = mibig.split('_')[0]
            if short_mibig in mibig_map:
                m = match(annotation, mibig_map[short_mibig])
                if m:
                    metadata = m
                    total_score += int(score)
                    print(m)
    return total_score, metadata

def aa_scoring(spectrum, gcf_like, tol=0.01):
    """
    Check for the prescence of AA mass shifts in the spectrum
    """
    from metabolomics import read_aa_losses
    aa_loss_file = 'aa_residues.csv'
    aa_losses = read_aa_losses(aa_loss_file)

    probs = []
    for bgc_aa_predictions in gcf_like.aa_predictions:
        p = 0
        for aa, aa_prob in bgc_aa_predictions.items():
            if aa_prob > 0:
                if aa in aa_losses:
                    mass_iso, mass_avg = aa_losses[aa]
                    found_losses = spectrum.has_loss(mass_iso, tol)
                    if len(found_losses) > 0:
                        p += aa_prob
        probs.append(p)

    return np.mean(probs)

def expand_spectrum_score(spectrum, gcf, scoring_function, strain_list):
    initial_score, initial_metadata = scoring_function(spectrum, gcf, strain_list)
    expanded_score, expanded_metadata = scoring_function(spectrum.family, gcf, strain_list)
    print("{} <-> {}\tInitial: {}, expanded: {} ({} spectra in family (id = {}))".format(
        spectrum, gcf,
        initial_score, expanded_score, len(spectrum.family.spectra),
        spectrum.family.family_id))
