import multiprocessing
import time
import itertools

import psutil
import numpy as np
import scipy.stats

# weights for metcalf_scoring, in order:
# - neither
# - GCF only
# - spectrum only
# - both
METCALF_WEIGHTS = [1, 0, -10, 10]

def compute_all_scores_multi(spectra_list, gcf_list, strain_list, scoring_function, do_random=True, cpus=8):
    spectra_part_list = np.array_split(spectra_list, cpus)

    # check actual number of available CPUs against requested number
    num_cpus = psutil.cpu_count()
    if cpus > num_cpus:
        print('Warning: {} CPUs requested, only {} available'.format(cpus, num_cpus))
        cpus = num_cpus

    q = multiprocessing.Queue()
    procs = []
    for i in range(cpus):
        if len(spectra_part_list[i]) == 0:
            continue

        p = multiprocessing.Process(target=compute_all_scores, args=(spectra_part_list[i], gcf_list, strain_list, scoring_function, do_random, q, i))
        procs.append(p)

    for p in procs:
        p.start()

    m_scores = {}

    t = time.time()
    num_finished = 0
    while num_finished < len(procs):
        data = q.get()
        m_scores.update(data)
        num_finished += 1

    print(('Total time: {:.1f}, {:.1f}/s'.format(time.time() - t, len(spectra_list) / (time.time() - t))))
    for p in procs:
        p.join()

    return m_scores

def compute_all_scores(spectra_list, gcf_list, strain_list, scoring_function, do_random=True, q=None, cpu_aff=None):
    m_scores = {}

    if cpu_aff is not None:
        psutil.Process().cpu_affinity([cpu_aff])

    for i, spectrum in enumerate(spectra_list):
        m_scores[spectrum] = {}
        if i % 100 == 0:
            print(("Done {} of {}".format(i, len(spectra_list))))
        for gcf in gcf_list:
            s, metadata = scoring_function(spectrum, gcf, strain_list)
            if do_random:
                s_random, _ = scoring_function(spectrum.random_spectrum, gcf.random_gcf, strain_list)
            else:
                s_random = None
            m_scores[spectrum][gcf] = (s, s_random, metadata)

    if q is not None:
        q.put(m_scores)

    return m_scores

def metcalf_scoring(spectral_like, gcf_like, strains, both=10, met_not_gcf=-10, gcf_not_met=0, neither=1):
    cum_score = 0
    shared_strains = set()
    for strain in strains:
        in_spec = spectral_like.has_strain(strain)
        in_gcf = gcf_like.has_strain(strain)
        if in_spec and in_gcf:
            cum_score += both
            shared_strains.add(strain)
        if in_spec and not in_gcf:
            cum_score += met_not_gcf
        if in_gcf and not in_spec:
            cum_score += gcf_not_met
        if not in_gcf and not in_spec:
            cum_score += neither
    return cum_score, shared_strains

def compute_all_scores_multi_np(spectra_list, gcf_list, strain_list, scoring_function, do_random=True, cpus=8):
    spectra_part_list = np.array_split(spectra_list, cpus)

    # check actual number of available CPUs against requested number
    num_cpus = psutil.cpu_count()
    if cpus > num_cpus:
        print('Warning: {} CPUs requested, only {} available'.format(cpus, num_cpus))
        cpus = num_cpus

    q = multiprocessing.Queue()

    procs = []
    for i in range(cpus):
        if len(spectra_part_list[i]) == 0:
            continue

        p = multiprocessing.Process(target=compute_all_scores_np, args=(spectra_part_list[i], gcf_list, strain_list, scoring_function, do_random, q, i))
        procs.append(p)

    for p in procs:
        p.start()

    t = time.time()
    num_finished = 0
    m_scores = {}

    # TODO it takes a *lot* of time to send all data back to the main process, almost 50% of 
    # the total when I was testing. Could probably be significantly reduced by avoiding
    # the use of Spectrum and GCF objects in m_scores to keep memory usage down and using
    # indices into an existing set of them held elsewhere.
    while num_finished < len(procs):
        data = q.get()
        m_scores.update(data)
        num_finished += 1

    for p in procs:
        p.join()

    print(('Total time: {:.1f} secs, {:.1f} scores/sec'.format(time.time() - t, len(spectra_list) / (time.time() - t))))
    return m_scores

def compute_all_scores_np(spectra_list, gcf_list, strain_list, scoring_function, do_random=True, q=None, cpu_aff=None):
    if cpu_aff is not None:
        psutil.Process().cpu_affinity([cpu_aff])

    if cpu_aff is not None:
        print('compute_all_scores_np on CPU {}, processing {} spectra'.format(cpu_aff, len(spectra_list)))
    else:
        print('compute_all_scores_np processing {} spectra'.format(len(spectra_list)))

    started = time.time()
    m_scores = {s: {} for s in spectra_list}

    # TODO: this could possibly be made faster by adjusting the datatypes slightly and using 
    # multiprocessing.Arrays instead, or a shared data structure via multiprocessing.Manager
    for spectrum, gcf in itertools.product(spectra_list, gcf_list):
        s, s_random, metadata = scoring_function(spectrum, gcf, strain_list, do_random)

        # don't care about results where there are no common strains, so exclude them here
        # (this saves a *lot* of entries in m_scores and makes the scoring significantly
        # faster overall)
        if metadata is not None:
            m_scores[spectrum][gcf] = (s, s_random, metadata)

    if q is not None:
        q.put(m_scores)

    ended = time.time()

    if cpu_aff is not None:
        print('compute_all_scores_np on CPU {}, total time = {:.1f}s, {:.1f} scores/sec'.format(cpu_aff, ended - started, len(m_scores) / (ended - started)))
    else:
        print('compute_all_scores_np: total time = {:.1f}s, {:.1f} scores/sec'.format(ended - started, len(m_scores) / (ended - started)))

    if cpu_aff is None:
        return m_scores

def metcalf_scoring_np(spectral_like, gcf_like, strains, do_random, weights=METCALF_WEIGHTS):
    # each strains_lookup array is treated as containing 1s/0s, so multiplying one
    # of them by 2 and adding together will give a new array containing values in
    # the range 0-3: 0 = neither, 1 = GCF only, 2 = spectrum only, 3 = both
    result = (2 * spectral_like.strains_lookup) + gcf_like.strains_lookup

    # use bincount to total up the number of occurrences of each value in result,
    # then just multiply by the weights and sum to get final score (minlength=4 ensures
    # that a 4-element array is always returned even if one of the 4 possible values
    # doesn't appear in result)
    cum_score = np.sum(np.bincount(result, minlength=4) * weights)

    cum_score_rand = 0
    if do_random:
        rand_result = (2 * spectral_like.random_spectrum.strains_lookup) + gcf_like.random_gcf.strains_lookup
        cum_score_rand = np.sum(np.bincount(rand_result, minlength=4) * weights)

    # get common strains, if any
    common = np.where(result == 3)[0]
    if len(common) > 0:
        return cum_score, cum_score_rand, strains[common]

    return cum_score, cum_score_rand, None

def hg_scoring(spectral_like, gcf_like, strains):
    spectral_count = 0
    gcf_count = 0
    overlap_count = 0

    for strain in strains:
        if spectral_like.has_strain(strain):
            spectral_count += 1
        if gcf_like.has_strain(strain):
            gcf_count += 1
        if spectral_like.has_strain(strain) and gcf_like.has_strain(strain):
            overlap_count += 1

    pos_in_sample = overlap_count
    N = spectral_count
    n = gcf_count
    M = len(strains)

    r = scipy.stats.hypergeom.sf(pos_in_sample, M, n, N, 1)
    return r, None

def name_scoring(spectral_like, gcf_like, mibig_map):
    score = 0
    metadata = None
    if len(spectral_like.annotations) == 0:
        print("No annotations")
        return None, None

    mibig_bgcs = gcf_like.get_mibig_bgcs()
    if len(mibig_bgcs) == 0:
        print("no mibig")
        return None, None

    for annotation in spectral_like.annotations:
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

def knownclusterblast_scoring(spectral_like, gcf_like, mibig_map):
    score = 0
    metadata = None
    if len(spectral_like.annotations) == 0:
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
    for annotation in spectral_like.annotations:
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
