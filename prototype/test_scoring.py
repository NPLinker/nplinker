import numpy as np

from nplinker.genomics import GCF, BGC
from nplinker.metabolomics import Spectrum, MolecularFamily, make_families
from nplinker.scoring.misc import metcalf_scoring, hg_scoring
from nplinker.data_linking import DataLinks, RandomisedDataLinks, LinkFinder
from nplinker.strains import Strain

np.random.seed(50)

def create_strains(n=10):
    return [Strain('strain_{:02d}'.format(x)) for x in range(n)]

def create_gcfs(strains, n=3):
    gcfs = []
    for i in range(n):
        gcf = GCF(i, 'fake_gcf_{}'.format(i), 'NRPS')
        num_strains = np.random.randint(1, len(strains))
        randoms = list(range(len(strains)))
        np.random.shuffle(randoms)
        randoms = randoms[:num_strains]
        gcf.id = i
        
        for j in range(num_strains):
            # (id, strain, name, bigscape_class, product_prediction, description=None):
            bgc = BGC(j, strains[randoms[j]], strains[randoms[j]].id, None, None)
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

def do_scoring_old(gcfs, spectra, strains):
    scores = {}
    for gcf in gcfs:
        scores[gcf] = {}
        for spec in spectra:
            score = metcalf_scoring(spec, gcf, strains)
            scores[gcf][spec] = score

    return scores

def do_scoring_new(gcfs, spectra, strains):
    datalinks = DataLinks()
    datalinks.load_data(spectra, gcfs, strains)
    datalinks.find_correlations()
    lf = LinkFinder()
    scores = lf.metcalf_scoring(datalinks)
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
    datalinks = DataLinks()
    datalinks.load_data(spectra, gcfs, strains)
    datalinks.find_correlations()
    lf = LinkFinder()
    scores = lf.hg_scoring(datalinks)
    return scores

def run_metcalf_test(n_strains=3, n_gcfs=5, n_spectra=4):
    strains = create_strains(n_strains)
    gcfs = create_gcfs(strains, n_gcfs)
    spectra = create_spectra(strains, n_spectra)

    old_scores = do_scoring_old(gcfs, spectra, strains)
    new_scores = do_scoring_new(gcfs, spectra, strains)

    print('Scores (old, new): ')
    for i, gcf in enumerate(gcfs):
        print('  {}: '.format(gcf))
        for j, spec in enumerate(spectra):
            print('    {}: {} | {} ({}/{})'.format(spec, old_scores[gcf][spec], new_scores[j][i], gcfs[i], spectra[j]))

if __name__ == "__main__":
    run_metcalf_test()
