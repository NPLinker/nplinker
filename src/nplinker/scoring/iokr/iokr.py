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

import numpy


class InputOutputKernelRegression:
    def __init__(self, data):
        self.data = data
        self.kernel_vector_cache = {}

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

    def calculate_fingerprint_kernel_vector(self, fingerprint, kernel="gaussian"):
        if kernel == "gaussian":

            def k(a, b, gamma=0.01):
                # kernel function
                d_sq = numpy.sum(numpy.power((a - b), 2))
                return numpy.exp(-gamma * d_sq)

            def k_vec(a_mat, b, gamma=0.01):
                # vectorised kernel function
                d_sq = numpy.sum(numpy.power((a_mat - b), 2), axis=1)
                return numpy.exp(-gamma * d_sq)

        # fp_id = hash(str(fingerprint))
        # if fp_id in self.kernel_vector_cache:
        #    return self.kernel_vector_cache[fp_id]

        training_data_latent = self.data.get_latent_vectors_vec(self.training_set)
        kernel_vector = k_vec(training_data_latent, fingerprint.T)
        # self.kernel_vector_cache[fp_id] = kernel_vector
        return kernel_vector

    def project_candidate(self, index, fingerprint):
        fingerprint_kernel_vector = self.calculate_fingerprint_kernel_vector(fingerprint)
        x_kernel_vector = self.data.kernel_product_set(index, self.training_set)
        res = numpy.dot(numpy.dot(fingerprint_kernel_vector, self.latent_basis), x_kernel_vector)
        return res

    def rank_candidates(self, index, candidate_fingerprints):
        candidate_distances = [
            self.project_candidate(index, fingerprint) for fingerprint in candidate_fingerprints
        ]
        return [
            x[1]
            for x in sorted(
                zip(candidate_distances, range(len(candidate_distances))),
                key=lambda x: x[0],
                reverse=True,
            )
        ]

    def project(self, index):
        x_kernel_vector = self.data.kernel_product_set(index, self.training_set)
        projection = numpy.dot(self.basis, x_kernel_vector)
        return projection

    def test(self, index, cutoff=0.01):
        proj = self.project(index)
        return [1 if x > cutoff else 0 for x in proj]


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
