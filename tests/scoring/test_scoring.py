import numpy as np
import pandas as pd
from nplinker.genomics import BGC
from nplinker.genomics import GCF
from nplinker.metabolomics.spectrum import Spectrum
from nplinker.scoring.linking.data_linking import DataLinks
from nplinker.scoring.linking.link_finder import LinkFinder
from nplinker.scoring.linking.misc_deprecated import hg_scoring
from nplinker.scoring.linking.misc_deprecated import metcalf_scoring
from nplinker.strains import Strain


np.random.seed(50)


def create_strains(n=10):
    return [Strain(f'strain_{x:02d}') for x in range(n)]


def create_gcfs(strains, n=3):
    gcfs = []
    for i in range(n):
        gcf = GCF(f'fake_gcf_{i}')
        num_strains = np.random.randint(1, len(strains))
        randoms = list(range(len(strains)))
        np.random.shuffle(randoms)
        randoms = randoms[:num_strains]

        for j in range(num_strains):
            bgc = BGC(j, strains[randoms[j]])
            gcf.add_bgc(bgc)

        gcfs.append(gcf)

    return gcfs


def create_spectra(strains, n=3):
    spectra = []

    families = np.random.randint(0, n * 2, n)

    for i in range(n):
        # (id, peaks, spectrum_id, precursor_mz, parent_mz=None, rt=None):
        spec = Spectrum(i, [(1, 2), (3, 4)], i, np.random.random())
        num_strains = np.random.randint(1, len(strains))
        randoms = list(range(len(strains)))
        np.random.shuffle(randoms)
        randoms = randoms[:num_strains]
        spec.family = int(families[i])
        for j in range(num_strains):
            spec.add_strain(strains[randoms[j]], 'foo', 1)

        spectra.append(spec)

    return spectra


def do_scoring_old(gcfs, spectra, strains, standardised):
    scores = {}
    for gcf in gcfs:
        scores[gcf] = {}
        for spec in spectra:
            score = metcalf_scoring(spec,
                                    gcf,
                                    strains,
                                    standardised=standardised)
            scores[gcf][spec] = score

    return scores


def do_scoring_new(gcfs, spectra, strains, standardised):
    datalinks = DataLinks(spectra, gcfs, strains)
    lf = LinkFinder()
    scores = lf.metcalf_scoring(datalinks)

    # print(lf.metcalf_expected)
    if standardised:
        # standardised scoring thing
        for i, gcf in enumerate(gcfs):
            for j, spec in enumerate(spectra):
                # get expected score, variance for objects with the current combo of strain counts
                # (note that spectrum = type 1 object here)
                met_strains = len(spec.strains)
                gen_strains = len(gcf.strains)
                expected = lf.metcalf_expected[met_strains][gen_strains]
                variance = lf.metcalf_variance[met_strains][gen_strains]
                scores[j][i] = (scores[j][i] - expected) / np.sqrt(variance)
    return scores


def do_scoring_old_hg(gcfs, spectra, strains):
    scores = {}
    for gcf in gcfs:
        scores[gcf] = {}
        for spec in spectra:
            score, _ = hg_scoring(spec, gcf, strains)
            scores[gcf][spec] = score

    return scores


def do_scoring_new_hg(gcfs, spectra, strains):
    datalinks = DataLinks(spectra, gcfs, strains)
    lf = LinkFinder()
    scores = lf.hg_scoring(datalinks)
    return scores


def run_metcalf_test(n_strains=3, n_gcfs=5, n_spectra=4, standardised=False):
    strains = create_strains(n_strains)
    gcfs = create_gcfs(strains, n_gcfs)
    spectra = create_spectra(strains, n_spectra)

    old_scores = do_scoring_old(gcfs, spectra, strains, standardised)
    new_scores = do_scoring_new(gcfs, spectra, strains, standardised)

    dfdata = {'nonvec_score': [], 'vec_score': [], 'gcf': [], 'spec': []}

    for i, gcf in enumerate(gcfs):
        for j, spec in enumerate(spectra):
            dfdata['nonvec_score'].append(old_scores[gcf][spec])
            dfdata['vec_score'].append(new_scores[j][i])
            dfdata['gcf'].append(gcf)
            dfdata['spec'].append(spec)

    return pd.DataFrame(data=dfdata)


if __name__ == "__main__":
    run_metcalf_test()
