import numpy
import os
import time

from . import iokr_opt
from . import iokrdata as iokrdataserver
from . import mk_fprints as fingerprint
from . import spectrum
from . import spectrum_filters

from logconfig import LogConfig
logger = LogConfig.getLogger(__file__)

class IOKRWrapper(object):
    """
    Wrapper around the IOKR server.
    Takes care of format conversion, fingerprint calculations, etc.
    Should also eventually take over the hardcoded stuff curently in get_iokr_server.
    """
    def __init__(self):
        self.fingerprint_type = None
        self.fingerprint_kernel = None
        self.ms_kernel = None

        self.iokr_server = None

    def _fingerprint(self, smiles):
        """
        Calculate molecular fingerprint for a SMILES string
        """
        return fingerprint.fingerprint_from_smiles(smiles, self.fingerprint_type)

    def rank_smiles(self, ms, candidate_smiles):
        """
        Rank a spectrum against a candidate set of SMILES strings
        """
        # TODO hacky
        spectrum_filters.datapath = get_datapath()
        ms.filter = spectrum_filters.filter_by_frozen_dag

        t = time.time()
        logger.debug('rank_smiles - Calculate candidate FPs')
        candidate_fps = numpy.array([self._fingerprint(x) for x in candidate_smiles])
        logger.debug('> {:.2f}s '.format(time.time() - t))
        t = time.time()
        logger.debug('rank_smiles - Extract latent basis')
        latent, latent_basis, gamma = self.iokr_server.get_data_for_novel_candidate_ranking()
        logger.debug('> {:.2f}s '.format(time.time() - t))
        t = time.time()
        logger.debug('rank_smiles - Get kernel vector for input sample')
        ms_kernel_vector = numpy.array(self.iokr_server.get_kernel_vector_for_sample(ms))
        logger.debug('> {:.2f}s '.format(time.time() - t))
        t = time.time()
        logger.debug('rank_smiles - Rank candidate set')
        ranking, _ = iokr_opt.rank_candidates_opt(0, candidate_fps, latent, ms_kernel_vector, latent_basis, gamma)
        logger.debug('> {:.2f}s '.format(time.time() - t))
        t = time.time()

        return ranking

def get_datapath():
    return os.path.join(os.path.dirname(__file__), 'data')

def get_iokr_server():
    datapath = get_datapath()
    kernel_files = [os.path.join(datapath, 'ppk_dag_all_normalised_shifted_{}.npy'.format(x)) for x in ['nloss', 'peaks']]

    fingerprint_type = "klekota-roth"
    iokr_wrapper = IOKRWrapper()

    iokr_wrapper.fingerprint_type = fingerprint_type
    iokr_wrapper.fingerprint_kernel = None  # function

    logger.debug('Init IOKR data server')
    iokrdata = iokrdataserver.IOKRDataServer(datapath, kernel=None)

    # When the kernel is initialised from matrix, we don't have guarantee
    # that the kernel_file and iokrdata.calculate_kernel match!
    logger.debug('Init kernel values')
    # Want to be able to set this to novel kernels
    kernel_matrix = iokrdataserver.load_kernels(kernel_files)
    iokrdata.kernel = kernel_matrix

    logger.debug('Load MS files')
    iokrdata.load_ms_files(datapath)

    def ppk_wrapper(ms_i, ms_j):
        sigma_mass = 0.00001
        sigma_int = 1000000.0
        ppk_peaks = spectrum.ppk(ms_i.spectrum, ms_j.spectrum, sigma_mass, sigma_int)
        ppk_nloss = spectrum.ppk_nloss(ms_i.spectrum, ms_j.spectrum, ms_i.parentmass, ms_j.parentmass, sigma_mass, sigma_int)

        if not hasattr(ms_i, 'ppk_peaks'):
            ms_i.ppk_peaks = spectrum.ppk(ms_i.spectrum, ms_i.spectrum, sigma_mass, sigma_int)
            ms_i.ppk_nloss = spectrum.ppk_nloss(ms_i.spectrum, ms_i.spectrum, ms_i.parentmass, ms_i.parentmass, sigma_mass, sigma_int)

        if not hasattr(ms_j, 'ppk_peaks'):
            ms_j.ppk_peaks = spectrum.ppk(ms_j.spectrum, ms_j.spectrum, sigma_mass, sigma_int)
            ms_j.ppk_nloss = spectrum.ppk_nloss(ms_j.spectrum, ms_j.spectrum, ms_j.parentmass, ms_j.parentmass, sigma_mass, sigma_int)

        ppk_peaks_normalised = ppk_peaks / numpy.sqrt(ms_i.ppk_peaks * ms_j.ppk_peaks)
        ppk_nloss_normalised = ppk_nloss / numpy.sqrt(ms_i.ppk_nloss * ms_j.ppk_nloss)
        ppk = ppk_peaks_normalised + ppk_nloss_normalised / 2
        return ppk

    # The function should accept two MSSpectrum objects and return a value
    logger.debug('Configure kernel')
    iokrdata.calculate_kernel = ppk_wrapper
    # TODO: This is super slow.
    # iokrdata.build_kernel_matrix()

    logger.debug('Set fingerprint')
    iokrdata.set_fingerprint(fingerprint_type)

    all_indices = iokrdata.get_all_indices()

    iokr = iokr_opt.InputOutputKernelRegression(iokrdata)
    logger.debug('fit()')
    iokr.set_training_indices(all_indices, _lambda=0.001)
    iokr.fit()

    iokr_wrapper.iokr_server = iokr

    logger.debug('get_iokr_server complete!')
    return iokr_wrapper

def test():
    from pyteomics import mgf
    # d = mgf.read(os.path.join(datapath, 'mibig/matched_mibig_gnps_2.0.mgf'))
    d = mgf.read(os.path.join(datapath, 'crusemann.mgf'))
    print(d)
    # Wrap the MGF entry in a MSSpectrum object
    test_spectrum = spectrum.MSSpectrum(d.next())

    # Candidate set
    SMILES = ["C\C=C\C=C\C(=O)C1=C(O)C(C)=C(O)C(C)=C1",
              "CC1=CC2=C(C(O)=C1)C(=O)C3=C(C=C(O)C=C3O)C2=O",
              "CC1=CC2=C(C(O)=C1)C(=O)C3=C(C=C(O)C=C3O)C2=O",
              "CC1=C2C(OC(=O)C3=C2C=C(O)C=C3O)=CC(O)=C1",
              "CC1=C2C(=O)C3=C(OC2=CC(O)=C1)C=C(O)C=C3O",
              "CC1CC(C)C(=O)C(C1)C(O)CC2CC(=O)NC(=O)C2"
              ]

    iokr = get_iokr_server()
    print('done init')

    print('rank')
    rank = iokr.rank_smiles(test_spectrum, SMILES)
    print(rank)


if __name__ == '__main__':
    test()
