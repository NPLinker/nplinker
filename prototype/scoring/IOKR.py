from metabolomics import Spectrum
from . iokr import nplinker_iokr
from . iokr.spectrum import MSSpectrum

from logconfig import LogConfig 
logger = LogConfig.getLogger(__file__)

def run_iokr_scoring(spec, bgc_list):
    logger.debug('IOKR datapath: {}'.format(nplinker_iokr.get_datapath()))

    print(spec)
    # extract SMILES strings for each supplied BGC (caching them in the object)
    smiles = [bgc.smiles for bgc in bgc_list]
    
    # TODO could probably handle this more gracefully by not submitting any BGC with no SMILES
    if None in smiles:
        logger.error('Failed to run IOKR scoring, one or more BGCs missing SMILES data')
        return None

    iokr_server = nplinker_iokr.get_iokr_server()
    rank = iokr_server.rank_smiles(MSSpectrum(spec=spec), smiles)
    return [bgc_list[i] for i in rank]

if __name__ == "__main__":
    # peaklist = [
	# (70.065002,84.632004),
	# (71.083000,39.488998),
	# (88.047997,54.071999),
	# (89.058998,471.227997),
	# (116.959999,75.615997),
	# (145.016998,68.410004),
	# (147.979996,62.558998),
	# (153.929993,94.302002),
	# (160.988007,261.783997),
	# (189.042999,389.872009),
	# (204.007996,43.602001),
	# (205.016006,159.628006),
    # ]
    # ms_test = Spectrum(id=0, peaks=peaklist, spectrum_id=0, precursor_mz=0, parent_mz=299.99701)
    run_iokr_scoring()
