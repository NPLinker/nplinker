# test functions

import numpy as np
from data_linking_functions import calc_correlation_matrix

def test_calc_correlation_matrix():
    # Test with simple input matrices and known outcome
    A = np.array([[0, 0, 0, 1], [0, 0, 1, 1], [0, 0, 1, 1], [0, 0, 0, 1], [0, 0, 0, 1]])
    B = np.array([[1, 1, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
    M1, M2, M3 = calc_correlation_matrix(A, B)
    assert M1[1][2] == 1    
    assert M2[2][1] == 2 
    assert M2[2][2] == 1 
    assert M3[4][0] == 2
    assert M3[1][2] == 0
    assert np.max([M1, M2, M3]) == 2  # in this example


from data_linking_functions import calc_likelihood_matrix
    
def test_calc_likelihood_matrix():
    # Test with simple input matrices and known outcome
    A = np.array([[0, 0, 0, 1], [0, 0, 1, 1], [0, 0, 1, 1], [0, 0, 0, 1], [0, 0, 0, 1]])
    B = np.array([[1, 1, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
    M1, M2, M3 = calc_correlation_matrix(A, B)
    LBA, LBnotA, LAB, LAnotB = calc_likelihood_matrix(A, B, M1, M2, M3)
    assert np.max([LBA, LBnotA, LAB, LAnotB]) <= 1  # no likelihood can be > 1
    assert LBA[1][2] == 0.5
    assert LAB[1][2] == 1
    assert LBnotA[1][2] == 0
    assert LBnotA[1][0] == 1
    assert LAnotB[1][0] == 1
    assert LAnotB[4][4] == 1/3
    assert LBA.shape == (len(A), len(B))  # must have shape len(A), len(B)
