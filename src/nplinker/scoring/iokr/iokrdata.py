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
import sys
import time
import numpy
import scipy.io
from . import mk_fprints
from . import spectrum
from . import spectrum_filters


logger = logging.getLogger(__name__)


def normalise_kernel(matrix):
    return matrix / numpy.sqrt(numpy.outer(matrix.diagonal(), matrix.diagonal()))


def load_kernel_file(filename, normalise=True):
    kernel = numpy.load(filename)
    if normalise:
        return normalise_kernel(kernel)
    else:
        return kernel


def load_kernels(kernel_files, normalise=True):
    kernel_matrices = [load_kernel_file(x, normalise) for x in kernel_files]
    kernel_sum = numpy.sum(kernel_matrices, axis=0)
    if normalise:
        return normalise_kernel(kernel_sum)
    else:
        return kernel_sum / len(kernel_matrices)


# Hold the GNPS records
class GNPS:
    def __init__(self, filename):
        logger.debug(f"GNPS {sys._getframe().f_code.co_name}")
        self.data_gnps = scipy.io.loadmat(filename)
        self.data_fp_array = numpy.array(self.data_gnps["fp"].todense())

    def get(self, index):
        # maybe add output type conversions (from sparse mat.)?
        fingerprint = self.data_fp_array[:, [index]]
        inchi = self.data_gnps["inchi"][index][0][0]
        formula = self.data_gnps["mf"][index][0][0]

        return inchi, formula, fingerprint

    def get_fingerprints(self, indices):
        logger.debug(f"GNPS {sys._getframe().f_code.co_name}")
        return self.data_fp_array[:, indices].T

    def set_fingerprint(self, fingerprint):
        logger.debug(f"GNPS {sys._getframe().f_code.co_name}")
        fp_list = []
        num_fps = len(self.data_gnps["inchi"])
        count = 0
        for inchi in self.data_gnps["inchi"]:
            count += 1
            inchi = inchi[0][0]
            fp = calc_fp(inchi, fingerprint)
            fp_list.append(fp)
            if count % 100 == 0:
                logger.debug(f"Done {count} / {num_fps}")
        self.data_fp_array = numpy.array(fp_list).T

    def set_fingerprint_from_file(self, filename):
        logger.debug(f"GNPS {sys._getframe().f_code.co_name}")
        self.data_fp_array = numpy.load(filename)

    def save_fingerprint_to_file(self, filename):
        numpy.save(filename, self.data_fp_array)


def load_folds(filename):
    with open(filename) as f:
        fold_ids = f.readlines()
    fold_ids = [x.strip() for x in fold_ids]
    return fold_ids


def load_kernel(filename):
    with open(filename) as f:
        raw = f.read().strip().split()

    data = numpy.array([float(x.strip()) for x in raw])
    dim = int(numpy.sqrt(len(data)))
    data = data.reshape((dim, dim))
    return data


def load_avg_kernel(path):
    kernel_files = [
        "ALIGND.txt",
        "CP2Plus.txt",
        "CPI.txt",
        "CPK.txt",
        "LB.txt",
        "LI.txt",
        "NI.txt",
        "RLB.txt",
        "ALIGN.txt",
        "CP2.txt",
        "CPJB.txt",
        "CSC.txt",
        "LC.txt",
        "LW.txt",
        "NSF.txt",
        "RLI.txt",
        "CEC.txt",
        "CPC.txt",
        "CPJ.txt",
        "FIPP.txt",
        "LIPP.txt",
        "NB.txt",
        "PPKr.txt",
        "WPC.txt",
    ]

    kernel = 0
    for i in kernel_files:
        kernel += load_kernel(path + os.sep + i)
    return kernel / len(kernel_files)


def load_candidates(path):
    candidate_sets = {}
    for file in os.listdir(path):
        candidate_name = extract_candidate_name(file)
        candidate_sets[candidate_name] = []
        if not file.endswith("mat"):
            continue
        data = scipy.io.loadmat(path + os.sep + file)
        num_candidates, _ = data["inchi"].shape
        for i in range(num_candidates):
            cand_inchi = data["inchi"][i, 0][0]
            cand_fingerprint = data["fp"][:, [i]]
            candidate_sets[candidate_name].append((cand_inchi, cand_fingerprint))

    return candidate_sets


def load_candidate_file_fp(filename, fp_filename, fingerprint):
    candidates = []
    data = scipy.io.loadmat(filename)
    num_candidates, _ = data["inchi"].shape
    cand_inchi = data["inchi"][:, 0]
    if os.path.exists(fp_filename + ".npy"):
        fp_data = numpy.load(fp_filename + ".npy")
    else:
        fp_data = [calc_fp(x[0], fingerprint) for x in cand_inchi]
        numpy.save(fp_filename, fp_data)
    for i in range(num_candidates):
        inchi = cand_inchi[i][0]
        cand_fingerprint = fp_data[i]
        candidates.append((inchi, cand_fingerprint))
    return candidates


def calc_fp(inchi, fingerprint):
    fp = mk_fprints.fingerprint_from_inchi(inchi, fingerprint)
    return fp


def load_candidate_file(filename):
    candidates = []
    data = scipy.io.loadmat(filename)
    num_candidates, _ = data["inchi"].shape
    # fp_vectors = numpy.array(data['fp'].todense())
    fp_vectors = numpy.array(data["fp"].todense())
    for i in range(num_candidates):
        cand_inchi = data["inchi"][i, 0][0]
        # cand_fingerprint = data['fp'][:, i]
        cand_fingerprint = fp_vectors[:, i]
        candidates.append((cand_inchi, cand_fingerprint))
    return candidates


def extract_candidate_name(filename):
    return filename.split(".")[0].split("_")[-1]


def load_spectra(filename):
    spectra = []
    with open(filename) as f:
        for l in f.readlines():
            if l.startswith("SPECTRUM_ID"):
                continue
            spectrum_id, compound, inchi = l.strip().split("\t")
            spectra.append((spectrum_id, compound, inchi))
    return spectra


# Holds the data that the IOKR can request.
# Needs to be central so we can get the kernel values for new samples
# TODO:
# - Also handle output kernel values (queriable in the same way?)
# - more flexible way of loading (input) kernels
# - Calclulate novel input kernels
class IOKRDataServer:
    def __init__(self, path, kernel=None):
        logger.debug("IOKRDataServer {}".format(sys._getframe().f_code.co_name))
        self.path = path
        self.gnps = GNPS(os.path.join(path, "data_GNPS.mat"))
        # self.candidates = load_candidates(path + os.sep + 'candidates')
        self.spectra = load_spectra(path + os.sep + "spectra.txt")

        self.folds = numpy.array(load_folds(path + os.sep + "cv_ind.txt"))
        if kernel is None:
            logger.debug("No kernel specified. Please initialise manually.")
            self.kernel = None
        elif kernel == "avg":
            logger.debug("Loading average kernel.")
            self.kernel = load_avg_kernel(path + os.sep + "input_kernels")
        else:
            logger.debug("Loading kernel %s" % kernel)
            self.kernel = load_kernel(path + os.sep + "input_kernels" + os.sep + kernel)

        self.candidate_path = path + os.sep + "candidates" + os.sep + "candidate_set_%s.mat"

        self.dimension = None
        self.fingerprint = None
        self.ms = []

    def load_ms_files(self, path):
        logger.debug("IOKRDataServer {}".format(sys._getframe().f_code.co_name))
        if os.path.exists(os.path.join(path, "ms.npy")):
            logger.debug("Loading cached MS data")
            t = time.time()
            try:
                self.ms = numpy.load(open(os.path.join(path, "ms.npy"), "rb"), allow_pickle=True)
                logger.debug(f"Done in {time.time() - t:2f}s")
                return
            except Exception as e:
                logger.debug("Failed to load cached MS data ({}), regenerating it...".format(e))

        logger.debug("Loading .ms files")
        t = time.time()
        for spectrum_id in [x[0] for x in self.spectra]:
            ms = spectrum.MSSpectrum()
            ms.correct_for_ionisation = True
            ms.normalise = True
            ms.filter = spectrum_filters.filter_by_frozen_dag
            ms.load(os.path.join(path, "SPEC", spectrum_id + ".ms"))
            self.ms.append(ms)
        logger.debug(f"Loading took {time.time() - t:2f}s")
        logger.debug("Caching MS data")
        numpy.save(open(os.path.join(path, "ms.npy"), "wb"), self.ms)

    def get_candidates(self, formula):
        logger.debug("IOKRDataServer {}".format(sys._getframe().f_code.co_name))
        if self.fingerprint is None:
            yield from load_candidate_file(self.candidate_path % formula)
        else:
            candidate_path = self.candidate_path % formula
            fp_path = self.fp_path + os.sep + "fp_%s.bin" % formula
            yield from load_candidate_file_fp(candidate_path, fp_path, self.fingerprint)

    def get_sample(self, idx, skip_candidates=False):
        gnps_inchi, formula, fingerprint = self.gnps.get(idx)
        if skip_candidates:
            candidates = []
        else:
            candidates = self.get_candidates(formula)
        spectrum_id, spectrum_name, spectrum_inchi = self.spectra[idx]

        assert gnps_inchi == spectrum_inchi
        return {
            "inchi": gnps_inchi,
            "formula": formula,
            "fingerprint": fingerprint,
            "candidates": candidates,
            "spectrum_id": spectrum_id,
            "spectrum_name": spectrum_name,
        }

    def get_cv_set(self, label, complement=False):
        logger.debug("IOKRDataServer {}".format(sys._getframe().f_code.co_name))
        if label not in self.folds:
            return []
        if not complement:
            indices = numpy.where(self.folds == label)[0]
        else:
            indices = numpy.where(self.folds != label)[0]
        kernel_submatrix = self.kernel[numpy.ix_(indices, indices)]
        sample = [self.get_sample(x, skip_candidates=True) for x in indices]
        fp_matrix = numpy.hstack([x["fingerprint"] for x in sample]).T

        return kernel_submatrix, fp_matrix

    def get_latent_vector(self, index):
        return (self.get_sample(index)["fingerprint"]).T

    def kernel_product(self, index_1, index_2):
        return self.kernel[index_1, index_2]

    def kernel_product_set(self, index_1, indices):
        return self.kernel[index_1, indices]

    def get_kernel_matrix(self, indices):
        return self.kernel[numpy.ix_(indices, indices)]

    def get_latent_vectors(self, indices):
        sample = [self.get_sample(x, skip_candidates=True) for x in indices]
        return numpy.hstack([x["fingerprint"] for x in sample]).T

    def get_latent_vectors_vec(self, indices):
        fingerprints = self.gnps.get_fingerprints(indices)
        return fingerprints

    def get_dimension(self):
        if self.dimension is None:
            self.dimension = len(self.get_sample(0, skip_candidates=True)["fingerprint"])
        return self.dimension

    def get_indices(self, label, complement=False):
        if complement:
            indices = numpy.where(self.folds != label)[0]
        else:
            indices = numpy.where(self.folds == label)[0]
        return indices

    def get_all_indices(self):
        return [x for x in range(len(self.folds))]

    def set_fingerprint(self, fingerprint):
        logger.debug("IOKRDataServer {}".format(sys._getframe().f_code.co_name))
        self.fingerprint = fingerprint

        fppath = self.path
        fpfile = os.path.join(self.path, f"fp_{fingerprint}_gnps.bin.npy")
        # fppath = self.path + os.sep + 'fp_' + fingerprint + os.sep
        # fpfile = fppath + os.sep + 'fp_gnps.bin'
        self.fp_path = fppath
        logger.debug(f'fpfile is "{fpfile}"')
        if os.path.exists(fpfile):
            logger.debug("Loading GNPS fingerprints from file")
            self.gnps.set_fingerprint_from_file(fpfile)
        else:
            logger.debug("Recalculating GNPS fingerprints")
            self.gnps.set_fingerprint(fingerprint)
            logger.debug("Saving GNPS fingerprints to file")
            self.gnps.save_fingerprint_to_file(fpfile)

    def build_kernel_matrix(self):
        logger.debug("IOKRDataServer {}".format(sys._getframe().f_code.co_name))
        kernel_matrix = numpy.zeros((len(self.spectra), len(self.spectra)))
        for i in range(len(self.spectra)):
            for j in range(i):
                kernel_value = self.calculate_kernel(self.ms[i], self.ms[j])
                kernel_matrix[i, j] = kernel_value
                kernel_matrix[j, i] = kernel_value

        self.kernel = normalise_kernel(kernel_matrix)
