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

import ctypes
import multiprocessing
from multiprocessing import Process
import numpy
from numba import jit


@jit(nopython=True)
def k_vec(a_mat, b, gamma=0.01):
    # vectorised kernel function
    d_sq = numpy.sum(numpy.power((a_mat - b), 2), axis=1)
    return numpy.exp(-gamma * d_sq)


@jit(nopython=True)
def project_candidate_opt(index, fingerprint, latent, x_kernel_vector, latent_basis, gamma):
    fingerprint_kernel_vector = k_vec(latent, fingerprint.T, gamma)
    res = numpy.dot(numpy.dot(fingerprint_kernel_vector, latent_basis), x_kernel_vector)
    return res, fingerprint_kernel_vector


@jit(nopython=True)
def preprocess_candidate_fingerprint(fingerprint, latent, latent_basis, gamma):
    fingerprint_kernel_vector = k_vec(latent, fingerprint.T, gamma)
    res = numpy.dot(fingerprint_kernel_vector, latent_basis)
    return res


@jit(nopython=True)
def project_candidate_preprocessed(preprocessed_fingerprint, x_kernel_vector):
    return numpy.dot(preprocessed_fingerprint, x_kernel_vector)


def project_candidate_wrapper(
    index, candidate_fingerprints, latent, x_kernel_vector, latent_basis, gamma
):
    return [
        project_candidate_opt(index, fingerprint, latent, x_kernel_vector, latent_basis, gamma)
        for fingerprint in candidate_fingerprints
    ]


@jit(nopython=True)
def project_candidate_kernel_opt(index, fingerprint_kernel_vector, x_kernel_vector, latent_basis):
    res = numpy.dot(numpy.dot(fingerprint_kernel_vector, latent_basis), x_kernel_vector)
    return res, fingerprint_kernel_vector


def project_candidate_kernel_wrapper(index, candidate_kernels, x_kernel_vector, latent_basis):
    return [
        project_candidate_kernel_opt(index, kernel_vector, x_kernel_vector, latent_basis)
        for kernel_vector in candidate_kernels
    ]


# @jit(nopython=True)
def rank_candidates_opt(
    index, candidate_fingerprints, latent, x_kernel_vector, latent_basis, gamma
):
    candidate_distances, fingerprint_kernel_vectors = project_candidates_opt(
        index, candidate_fingerprints, latent, x_kernel_vector, latent_basis, gamma
    )
    return numpy.argsort(numpy.array(candidate_distances))[::-1], fingerprint_kernel_vectors


def project_candidates_opt(
    index, candidate_fingerprints, latent, x_kernel_vector, latent_basis, gamma
):
    candidate_projections = [
        project_candidate_opt(index, fingerprint, latent, x_kernel_vector, latent_basis, gamma)
        for fingerprint in candidate_fingerprints
    ]
    candidate_distances = [x[0] for x in candidate_projections]
    fingerprint_kernel_vectors = [x[1] for x in candidate_projections]
    return candidate_distances, fingerprint_kernel_vectors


def preprocess_candidates(candidate_fingerprints, latent, latent_basis, gamma):
    fingerprint_projections = [
        preprocess_candidate_fingerprint(fp, latent, latent_basis, gamma)
        for fp in candidate_fingerprints
    ]
    return fingerprint_projections


def project_candidates_preprocessed(preprocessed_candidates, x_kernel_vector):
    candidate_distances = [
        project_candidate_preprocessed(x, x_kernel_vector) for x in preprocessed_candidates
    ]
    return candidate_distances, []


def rank_candidate_kernel_opt(index, candidate_kernels, latent, x_kernel_vector, latent_basis):
    candidate_projections = [
        project_candidate_kernel_opt(index, kernel_vector, x_kernel_vector, latent_basis)
        for kernel_vector in candidate_kernels
    ]
    candidate_distances = [x[0] for x in candidate_projections]
    fingerprint_kernel_vectors = [x[1] for x in candidate_projections]
    return numpy.argsort(numpy.array(candidate_distances))[::-1], fingerprint_kernel_vectors


class InputOutputKernelRegression:
    def __init__(self, data):
        self.data = data
        self.kernel_vector_cache = {}

        self.ker_v_list = []

    def set_training_indices(self, indices=None, _lambda=0.1):
        self._lambda = _lambda
        if indices is None:
            indices = range(self.data.data_size)
        self.training_set = indices

    def fit(self):
        training_data_kernel = self.data.get_kernel_matrix(self.training_set)
        training_data_latent = self.data.get_latent_vectors(self.training_set)

        eye = numpy.eye(len(training_data_kernel))

        training_data_latent = numpy.array(training_data_latent).T

        latent_basis = numpy.linalg.inv(self._lambda * eye + training_data_kernel)
        self.latent_basis = latent_basis
        self.basis = numpy.dot(training_data_latent, latent_basis)

        self.latent = self.data.get_latent_vectors_vec(self.training_set)

    def calculate_fingerprint_kernel_vector(
        self, fingerprint, training_data_latent, kernel="gaussian"
    ):
        # if kernel == 'gaussian':
        #     def k_vec(a_mat, b, gamma=0.01):
        #         # vectorised kernel function
        #         d_sq = numpy.sum(numpy.power((a_mat - b), 2), axis=1)
        #         return numpy.exp(- gamma * d_sq)

        kernel_vector = k_vec(training_data_latent, fingerprint.T)
        return kernel_vector

    def project_candidate(self, index, fingerprint, latent, output=None):
        fingerprint_kernel_vector = self.calculate_fingerprint_kernel_vector(fingerprint, latent)
        x_kernel_vector = self.data.kernel_product_set(index, self.training_set)
        res = numpy.dot(numpy.dot(fingerprint_kernel_vector, self.latent_basis), x_kernel_vector)

        if output is not None:
            output.value = res
        return res

    def rank_candidates(self, index, candidate_fingerprints):
        latent = self.data.get_latent_vectors_vec(self.training_set)

        # candidate_distances = [self.project_candidate(index, fingerprint, latent) for fingerprint in candidate_fingerprints]

        candidate_distances = []
        batch_size = 100
        for i in range(int(len(candidate_fingerprints) / batch_size) + 1):
            if i % 10 == 0:
                print("processing set %s" % i)
            candidate_fingerprints_subset = candidate_fingerprints[
                i * batch_size : (i + 1) * batch_size
            ]
            values = [
                multiprocessing.Value(ctypes.c_float)
                for x in range(len(candidate_fingerprints_subset))
            ]
            candidate_distance_jobs = [
                Process(target=self.project_candidate, args=(index, fingerprint, latent, value))
                for fingerprint, value in zip(candidate_fingerprints_subset, values)
            ]
            [x.start() for x in candidate_distance_jobs]
            [x.join() for x in candidate_distance_jobs]
            candidate_distances_subset = [x.value for x in values]
            candidate_distances.extend(candidate_distances_subset)

        return [
            x[1]
            for x in sorted(
                zip(candidate_distances, range(len(candidate_distances))),
                key=lambda x: x[0],
                reverse=True,
            )
        ]

    def rank_candidates_opt(self, index, candidate_fingerprints, candidate_set_name=None):
        latent = self.latent
        x_kernel_vector = self.data.kernel_product_set(index, self.training_set)
        latent_basis = self.latent_basis
        gamma = 0.01

        candidate_ranking, fp_kernels = rank_candidates_opt(
            index, candidate_fingerprints, latent, x_kernel_vector, latent_basis, gamma
        )

        # if candidate_set_name is not None:
        #     with open('/home/grimur/iokr/data/candidate_hash/fp_gamma0.01_%s.bin' % candidate_set_name, 'wb') as f:
        #         pickle.dump(fp_kernels, f)

        return candidate_ranking

    def get_data_for_candidate_ranking(self, index):
        latent, latent_basis, gamma = self.get_data_for_novel_candidate_ranking()
        x_kernel_vector = self.data.kernel_product_set(index, self.training_set)

        return latent, x_kernel_vector, latent_basis, gamma

    def get_data_for_novel_candidate_ranking(self):
        latent = self.latent
        latent_basis = self.latent_basis
        gamma = 0.01
        return latent, latent_basis, gamma

    def project(self, index):
        x_kernel_vector = self.data.kernel_product_set(index, self.training_set)
        projection = numpy.dot(self.basis, x_kernel_vector)
        return projection

    def test(self, index, cutoff=0.01):
        proj = self.project(index)
        return [1 if x > cutoff else 0 for x in proj]

    def get_kernel_vector_for_sample(self, ms):
        kernel_vector = []

        ms_auto = self.data.calculate_kernel(ms, ms)
        for idx in self.training_set:
            t_sp = self.data.ms[idx]
            t_sp_auto = self.data.calculate_kernel(t_sp, t_sp)
            ms_t_sp = self.data.calculate_kernel(ms, t_sp)
            kernel_value = ms_t_sp / numpy.sqrt(ms_auto * t_sp_auto)
            kernel_vector.append(kernel_value)

        return kernel_vector


def main():
    target_vectors = numpy.array(
        [[0, 1, 0], [1, 0, 0], [1, 1, 0], [1, 1, 1], [0, 0, 1], [0, 1, 1], [1, 0, 1]], dtype="float"
    )
    repr_vectors = numpy.array(
        [[0, 1, 0], [1, 0, 0], [1, 1, 0], [1, 1, 1], [0, 0, 1], [0, 1, 1], [1, 0, 1]], dtype="float"
    )

    kernel_matrix = numpy.zeros((7, 7))
    for i in range(7):
        for j in range(i + 1):
            kernel_matrix[i, j] = kernel_matrix[j, i] = numpy.dot(repr_vectors[i], repr_vectors[j])

    data = DataStore(kernel_matrix, target_vectors)
    okr = InputOutputKernelRegression(data)
    okr.set_training_indices([0, 2, 3], _lambda=0.0)
    okr.fit()
    print(okr.test(0), target_vectors[0])
    print(okr.test(1), target_vectors[1])
    print(okr.test(2), target_vectors[2])
    print(okr.test(3), target_vectors[3])


if __name__ == "__main__":
    main()
