# Copyright 2021 The NPLinker Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import logging
import os
import time
import numpy
from nplinker.genomics import GCF
from nplinker.metabolomics import MolecularFamily
from . import iokr_opt
from . import iokrdata as iokrdataserver
from . import mk_fprints as fingerprint
from . import spectrum
from . import spectrum_filters


logger = logging.getLogger(__name__)


class NPLinkerIOKR:
    """Rank the potential links in a NPLinker object using IOKR, and wrap the results."""

    def __init__(self, npl):
        self.npl = npl

        # Build BGC<->SMILES lookup
        bgc_smiles_lookup = {}
        for bgc_idx, bgc in enumerate(npl.bgcs):
            bgc_smiles_lookup[bgc] = []
            if bgc.smiles is not None:
                try:
                    fingerprint.fingerprint_from_smiles(bgc.smiles)
                    bgc_smiles_lookup[bgc].append(bgc.smiles)
                except ValueError:
                    print(f"Filtered out smiles {bgc.smiles}")
        self.bgc_smiles_lookup = bgc_smiles_lookup

        smiles_bgc_lookup = {}
        for bgc, smiles in bgc_smiles_lookup.items():
            if smiles is None:
                continue
            for smiles_string in smiles:
                if smiles_string in smiles_bgc_lookup:
                    smiles_bgc_lookup[smiles_string].add(bgc)
                else:
                    smiles_bgc_lookup[smiles_string] = {bgc}

        self.smiles_list = list(smiles_bgc_lookup.keys())

        # Map NPLinker spectra to format understood by IOKR
        self.spectra_list = [spectrum.MSSpectrum(spec=x) for x in npl.spectra]

        # Initialise IOKR server
        self.iokr_server = get_iokr_server()

    def calculate_scores(self):
        self.iokr_scores = self.iokr_server.score_smiles(self.spectra_list, self.smiles_list)

    def get_score(self, mf_like, gcf_like):
        if isinstance(mf_like, MolecularFamily):
            mf = mf_like
        else:
            mf = MolecularFamily(-1)
            mf.spectra = [mf_like]

        if isinstance(gcf_like, GCF):
            gcf = gcf_like
        else:
            gcf = GCF("-1")
            gcf.bgcs = [gcf_like]

        potential_scores = []
        for spec in mf.spectra:
            spectrum_idx = self.spectra_list.index(spec)
            for bgc in gcf.bgcs:
                bgc_smiles_list = self.bgc_smiles_lookup[bgc]
                for smiles in bgc_smiles_list:
                    smiles_idx = self.smiles_list.index(smiles)
                    spectrum_smiles_score = self.iokr_scores[spectrum_idx, smiles_idx]
                    potential_scores.append(spectrum_smiles_score)

        if len(potential_scores) == 0:
            logger.info(f"No scores found for {gcf_like}")
            return 0
        else:
            return max(potential_scores)


class IOKRWrapper:
    """Wrapper around the IOKR server.
    Takes care of format conversion, fingerprint calculations, etc.
    Should also eventually take over the hardcoded stuff curently in get_iokr_server.
    """

    def __init__(self):
        self.fingerprint_type = None
        self.fingerprint_kernel = None
        self.ms_kernel = None

        self.iokr_server = None

    def _fingerprint(self, smiles):
        """Calculate molecular fingerprint for a SMILES string."""
        return fingerprint.fingerprint_from_smiles(smiles, self.fingerprint_type)

    def score_smiles(self, ms_list, candidate_smiles):
        """Score a set of spectra against a candidate set of SMILES strings."""
        spectrum_filters.datapath = get_datapath()

        logger.debug("cache miss")
        logger.debug("Calculate fingerprints for candidate set")
        # TODO: Cache this
        candidate_fps = []
        for i, c in enumerate(candidate_smiles):
            logger.debug(f"done {i}/{len(candidate_smiles)}")
            candidate_fps.append(self._fingerprint(c))
        candidate_fps = numpy.array(candidate_fps)
        logger.debug("Extract latent basis")
        latent, latent_basis, gamma = self.iokr_server.get_data_for_novel_candidate_ranking()
        logger.debug("writing cache")

        projection_matrix = numpy.zeros((len(ms_list), len(candidate_smiles)))

        # TODO: Cache this
        logger.debug("Preprocessing candidate set FPs")
        candidates = iokr_opt.preprocess_candidates(candidate_fps, latent, latent_basis, gamma)

        for ms_index, ms in enumerate(ms_list):
            logger.debug("Rank spectrum {} ({}/{})".format(ms.id, ms_index, len(ms_list)))
            ms.filter = spectrum_filters.filter_by_frozen_dag
            logger.debug("kernel vector")
            t0 = time.time()
            ms_kernel_vector = numpy.array(self.iokr_server.get_kernel_vector_for_sample(ms))
            t1 = time.time()
            logger.debug(f"done ({t1 - t0})")
            logger.debug("project")
            # projections, _ = iokr_opt.project_candidates_opt(0, candidate_fps, latent, ms_kernel_vector, latent_basis, gamma)
            projections, _ = iokr_opt.project_candidates_preprocessed(candidates, ms_kernel_vector)
            t2 = time.time()
            logger.debug(f"done ({t2 - t1})")

            logger.debug("save distances")
            projection_matrix[ms_index, :] = projections

        return projection_matrix

    def rank_smiles(self, ms, candidate_smiles):
        """Rank a spectrum against a candidate set of SMILES strings."""
        # TODO hacky
        spectrum_filters.datapath = get_datapath()
        ms.filter = spectrum_filters.filter_by_frozen_dag

        t = time.time()
        logger.debug("rank_smiles - Calculate candidate FPs")
        candidate_fps = numpy.array([self._fingerprint(x) for x in candidate_smiles])
        logger.debug(f"> {time.time() - t:.2f}s ")
        t = time.time()
        logger.debug("rank_smiles - Extract latent basis")
        latent, latent_basis, gamma = self.iokr_server.get_data_for_novel_candidate_ranking()
        logger.debug(f"> {time.time() - t:.2f}s ")
        t = time.time()
        logger.debug("rank_smiles - Get kernel vector for input sample")
        ms_kernel_vector = numpy.array(self.iokr_server.get_kernel_vector_for_sample(ms))
        logger.debug(f"> {time.time() - t:.2f}s ")
        t = time.time()
        logger.debug("rank_smiles - Rank candidate set")
        ranking, _ = iokr_opt.rank_candidates_opt(
            0, candidate_fps, latent, ms_kernel_vector, latent_basis, gamma
        )
        logger.debug(f"> {time.time() - t:.2f}s ")
        t = time.time()

        return ranking


def get_datapath():
    return os.path.join(os.path.dirname(__file__), "data")


def get_iokr_server():
    datapath = get_datapath()
    kernel_files = [
        os.path.join(datapath, f"ppk_dag_all_normalised_shifted_{x}.npy")
        for x in ["nloss", "peaks"]
    ]

    # fingerprint_type = "klekota-roth"
    fingerprint_type = None
    iokr_wrapper = IOKRWrapper()

    iokr_wrapper.fingerprint_type = fingerprint_type
    iokr_wrapper.fingerprint_kernel = None  # function

    logger.debug("Init IOKR data server")
    iokrdata = iokrdataserver.IOKRDataServer(datapath, kernel=None)

    # When the kernel is initialised from matrix, we don't have guarantee
    # that the kernel_file and iokrdata.calculate_kernel match!
    logger.debug("Init kernel values")
    # Want to be able to set this to novel kernels
    kernel_matrix = iokrdataserver.load_kernels(kernel_files)
    iokrdata.kernel = kernel_matrix

    logger.debug("Load MS files")
    iokrdata.load_ms_files(datapath)

    def ppk_wrapper(ms_i, ms_j):
        sigma_mass = 0.00001
        sigma_int = 1000000.0
        ppk_peaks = spectrum.ppk(ms_i.spectrum, ms_j.spectrum, sigma_mass, sigma_int)
        ppk_nloss = spectrum.ppk_nloss(
            ms_i.spectrum, ms_j.spectrum, ms_i.parentmass, ms_j.parentmass, sigma_mass, sigma_int
        )

        if not hasattr(ms_i, "ppk_peaks"):
            ms_i.ppk_peaks = spectrum.ppk(ms_i.spectrum, ms_i.spectrum, sigma_mass, sigma_int)
            ms_i.ppk_nloss = spectrum.ppk_nloss(
                ms_i.spectrum,
                ms_i.spectrum,
                ms_i.parentmass,
                ms_i.parentmass,
                sigma_mass,
                sigma_int,
            )

        if not hasattr(ms_j, "ppk_peaks"):
            ms_j.ppk_peaks = spectrum.ppk(ms_j.spectrum, ms_j.spectrum, sigma_mass, sigma_int)
            ms_j.ppk_nloss = spectrum.ppk_nloss(
                ms_j.spectrum,
                ms_j.spectrum,
                ms_j.parentmass,
                ms_j.parentmass,
                sigma_mass,
                sigma_int,
            )

        ppk_peaks_normalised = ppk_peaks / numpy.sqrt(ms_i.ppk_peaks * ms_j.ppk_peaks)
        ppk_nloss_normalised = ppk_nloss / numpy.sqrt(ms_i.ppk_nloss * ms_j.ppk_nloss)
        ppk = ppk_peaks_normalised + ppk_nloss_normalised / 2
        return ppk

    # The function should accept two MSSpectrum objects and return a value
    logger.debug("Configure kernel")
    iokrdata.calculate_kernel = ppk_wrapper
    # TODO: This is super slow.
    # iokrdata.build_kernel_matrix()

    logger.debug("Set fingerprint")
    iokrdata.set_fingerprint(fingerprint_type)

    all_indices = iokrdata.get_all_indices()

    iokr = iokr_opt.InputOutputKernelRegression(iokrdata)
    logger.debug("fit()")
    iokr.set_training_indices(all_indices, _lambda=0.001)
    iokr.fit()

    iokr_wrapper.iokr_server = iokr

    logger.debug("get_iokr_server complete!")
    return iokr_wrapper


def test():
    from pyteomics import mgf

    # d = mgf.read(os.path.join(datapath, 'mibig/matched_mibig_gnps_2.0.mgf'))
    d = mgf.read(os.path.join(datapath, "crusemann.mgf"))
    print(d)
    # Wrap the MGF entry in a MSSpectrum object
    test_spectrum = spectrum.MSSpectrum(d.next())

    # Candidate set
    SMILES = [
        r"C\C=C\C=C\C(=O)C1=C(O)C(C)=C(O)C(C)=C1",
        "CC1=CC2=C(C(O)=C1)C(=O)C3=C(C=C(O)C=C3O)C2=O",
        "CC1=CC2=C(C(O)=C1)C(=O)C3=C(C=C(O)C=C3O)C2=O",
        "CC1=C2C(OC(=O)C3=C2C=C(O)C=C3O)=CC(O)=C1",
        "CC1=C2C(=O)C3=C(OC2=CC(O)=C1)C=C(O)C=C3O",
        "CC1CC(C)C(=O)C(C1)C(O)CC2CC(=O)NC(=O)C2",
    ]

    iokr = get_iokr_server()
    print("done init")

    print("rank")
    rank = iokr.rank_smiles(test_spectrum, SMILES)
    print(rank)


if __name__ == "__main__":
    test()
