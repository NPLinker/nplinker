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

# test functions

import numpy as np
from nplinker.scoring.linking.data_linking_functions import calc_correlation_matrix


def test_calc_correlation_matrix():
    # Test with simple input matrices and known outcome
    A = np.array([[0, 0, 0, 1], [0, 0, 1, 1], [0, 0, 1, 1], [0, 0, 0, 1],
                  [0, 0, 0, 1]])
    B = np.array([[1, 1, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0],
                  [0, 1, 0, 0], [0, 0, 0, 1]])
    M_A_B, M_A_notB, M_notA_B, M_notA_notB = calc_correlation_matrix(A, B)
    assert M_A_B[1][2] == 1
    assert M_A_notB[2][1] == 2
    assert M_A_notB[2][2] == 1
    assert M_notA_B[4][0] == 2
    assert M_notA_B[1][2] == 0
    assert np.max([M_A_B, M_A_notB, M_notA_B]) == 2  # in this example
    assert np.max(M_notA_notB) == 3
    assert np.max(M_notA_notB == np.array(
        [[1, 2, 2, 2, 2, 3], [0, 1, 2, 1, 1, 2], [0, 1, 2, 1, 1, 2],
         [1, 2, 2, 2, 2, 3], [1, 2, 2, 2, 2, 3]]))


from nplinker.scoring.linking.data_linking_functions import calc_likelihood_matrix


def test_calc_likelihood_matrix():
    # Test with simple input matrices and known outcome
    A = np.array([[0, 0, 0, 1], [0, 0, 1, 1], [0, 0, 1, 1], [0, 0, 0, 1],
                  [0, 0, 0, 1]])
    B = np.array([[1, 1, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0],
                  [0, 1, 0, 0], [0, 0, 0, 1]])
    M_A_B, M_A_notB, M_notA_B, M_notA_notB = calc_correlation_matrix(A, B)
    LBA, LBnotA, LAB, LAnotB = calc_likelihood_matrix(A, B, M_A_B, M_A_notB,
                                                      M_notA_B)
    assert np.max([LBA, LBnotA, LAB, LAnotB]) <= 1  # no likelihood can be > 1
    assert LBA[1][2] == 0.5
    assert LAB[1][2] == 1
    assert LBnotA[1][2] == 0
    assert LBnotA[1][0] == 1
    assert LAnotB[1][0] == 1
    assert LAnotB[4][4] == 1 / 3
    assert LBA.shape == (len(A), len(B))  # must have shape len(A), len(B)


from nplinker.scoring.linking.data_linking_functions import pair_prob_hg


def test_pair_prob_hg():
    # Test pair_prob with known cases
    assert pair_prob_hg(1, 100, 1, 1) == 1 / 100
    assert pair_prob_hg(1, 100, 50, 1) == 0.5
    assert pair_prob_hg(1, 100, 1, 50) == 0.5
    assert pair_prob_hg(1, 100, 2, 2) == 98 / 100 * 2 / 99 + 2 / 100 * 98 / 99


from nplinker.scoring.linking.data_linking_functions import hit_prob_dist


def test_hit_prob_dist():
    # Testhit_prob_dist with known cases

    pks = hit_prob_dist(100, 1, 1, 100)
    assert np.sum(pks) > 0.99999999
    assert np.sum(pks) < 1.00000001
    assert pks[0][0] == 0.99**100


from nplinker.scoring.linking.data_linking_functions import permutation_unique


def test_permutation_unique():
    # Test function to find unique permutations with known cases

    testlist = [1, 2, 3, 4, 5]
    assert len(list(permutation_unique(testlist))) == 120  # math.factorial(5)
    testlist = [1, 2, 3, 4, 1]
    assert len(list(permutation_unique(
        testlist))) == 60  # math.factorial(5)/math.factorial(2)
    testlist = [1, 2, 3, 1, 2, 3]
    assert len(list(permutation_unique(
        testlist))) == 90  # math.factorial(6)/math.factorial(2)**3
    testlist = ['A', 'B', 'C', 'C', 'C']
    assert len(list(permutation_unique(
        testlist))) == 20  # math.factorial(5)/math.factorial(3)


from nplinker.scoring.linking.data_linking_functions import pair_prob


def test_pair_prob():
    # Test pair_prob function with known cases
    P_str = np.ones(10)
    P_str = P_str / np.sum(P_str)
    XG = [0, 1, 2]
    Ny = 3
    hits = 2
    assert pair_prob(P_str, XG, Ny, hits) == 0.175

    P_str = np.ones(10)
    P_str = P_str / np.sum(P_str)
    XG = [0, 1]
    Ny = 2
    hits = 2
    assert pair_prob(
        P_str, XG, Ny, hits
    ) > 2 / 90  # correct result here actually should be = 2/90, but there are rounding errors
    assert pair_prob(P_str, XG, Ny, hits) < (2 / 90 + 0.00000001)


from nplinker.scoring.linking.data_linking_functions import link_prob


def test_link_prob():
    # Test link_prob function with known cases
    Nstr = 10
    P_str = np.ones(Nstr)
    P_str[0:2] = P_str[0:2] * 0.1
    P_str = P_str / np.sum(P_str)
    XGS = [0, 1]
    Ny = 2
    Nx = 2
    assert link_prob(P_str, XGS, Nx, Ny, Nstr) == 2 * P_str[1]**2 / 0.9

    Ny = 3
    Nx = 5
    assert link_prob(P_str, XGS, Nx, Ny,
                     Nstr) == 6 * 0.5 * P_str[0]**2 / (0.9 * 0.8)
