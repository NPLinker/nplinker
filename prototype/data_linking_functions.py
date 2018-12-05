# functions

import numpy as np

def calc_correlation_matrix(M_type1_cond, M_type2_cond):
    """ 
    Calculate correlation matrices from co-occurence matrices
    Input:
    M_type1_cond(x,y) is 1 if type1_x IS observed under condition_y
    M_type1_cond(x,y) is 0 if type1_x IS NOT observed under condition_y
    
    Outputs three correlation matrices:
    M_type1_type2(x,y) --- number of conditions where type1_x and type2_y co-occur
    M_type1_nottype2(x,y) --- number of conditions where type1_x and NOT-type2_y co-occur
    M_nottype1_type2(x,y) --- number of conditions where NOT-type1_x and type2_y co-occur 
    """
          
    # Quick computation of sum both present
    testA = np.dot(M_type1_cond,M_type2_cond.T)
    # Present in type1 and not in type 2
    testB = np.dot(M_type1_cond,1-M_type2_cond.T)
    # Not in type 1 and in type 2
    testC = np.dot(1-M_type1_cond,M_type2_cond.T)
    # Not in either
    testD = np.dot(1-M_type1_cond,1-M_type2_cond.T)
    # print(np.abs((testA - M_type1_type2)).sum())
    # print(np.abs((testB - M_type1_nottype2)).sum())
    # print(np.abs((testC - M_nottype1_type2)).sum())
    M_type1_type2 = testA
    M_type1_nottype2 = testB
    M_nottype1_type2 = testC
    M_nottype1_nottype2 = testD
    return M_type1_type2, M_type1_nottype2, M_nottype1_type2 , M_nottype1_nottype2



def calc_likelihood_matrix(M_type1_cond, M_type2_cond, 
                           M_type1_type2, M_type1_nottype2, M_nottype1_type2):
    """ 
    Calculate correlation matrices from co-occurence matrices
    Input:
    M_type1_cond(x,y) is 1 if type1_x IS observed under condition_y
    M_type1_cond(x,y) is 0 if type1_x IS NOT observed under condition_y
    M_type1_type2(x,y) --- number of conditions where type1_x and type2_y co-occur
    M_type1_nottype2(x,y) --- number of conditions where type1_x and NOT-type2_y co-occur
    M_nottype1_type2(x,y) --- number of conditions where NOT-type1_x and type2_y co-occur

    Output:
    Four likelihood matrices of size len(type1) x len(type2):
    P_type2_given_type1, P_type2_not_type1, P_type1_given_type2, P_type1_not_type2
    """
    
    dim1, dim2 = M_type1_type2.shape
    num_conditions = M_type2_cond.shape[1]

    P_type2_given_type1 = np.zeros((dim1, dim2))
    P_type2_not_type1 = np.zeros((dim1, dim2))
    P_type1_given_type2 = np.zeros((dim1, dim2))
    P_type1_not_type2 = np.zeros((dim1, dim2))

    # Calculate P_type2_given_type1 matrix
    P_sum_type1 = np.sum(M_type1_cond, axis=1)
    P_sum_type1[P_sum_type1 < 1] = 1 #avoid later division by 0
    P_type2_given_type1 = M_type1_type2/np.tile(P_sum_type1, (dim2, 1)).transpose(1, 0)

    # Calculate P_type1_given_type2 matrix
    P_sum_type2 = np.sum(M_type2_cond, axis=1)
    P_sum_type2[P_sum_type2 < 1] = 1 #avoid later division by 0
    P_type1_given_type2 = M_type1_type2/np.tile(P_sum_type2, (dim1, 1))

    # Calculate P_type2_not_type1 matrix
    P_sum_not_type1 = num_conditions - P_sum_type1
    P_type2_not_type1 = M_nottype1_type2/np.tile(P_sum_not_type1, (dim2, 1)).transpose(1, 0)

    # Calculate P_type1_not_type2 matrix
    P_sum_not_type2 = num_conditions - P_sum_type2  
    P_type1_not_type2 = M_type1_nottype2/np.tile(P_sum_not_type2, (dim1, 1))
    
    return P_type2_given_type1, P_type2_not_type1, P_type1_given_type2, P_type1_not_type2


def pair_prob(k, N, Nx, Ny):
    """
    Calculate the probability to draw k times type(Ny) out of N elements (whereof Ny type(Ny)s),
    when drawing Nx times in total.
    """
    if (k > Nx) or (k > Ny):
        print("Given 'k' must be <= Nx and <= Ny.")
    
    import math
    
    term1 = math.factorial(Ny)/(math.factorial(k) * math.factorial(Ny - k))
    term2 = math.factorial(N - Ny)/(math.factorial(Nx - k) * math.factorial(N - Nx - Ny + k))
    term3 = math.factorial(N)/(math.factorial(Nx) * math.factorial(N - Nx))
    
    return term1 * term2 / term3 
    

def hit_prob_dist(N, Nx, Ny, nys):
    
    import math
    
    p_dist_ks = []
    for k in range(1,min(Nx, Ny)+1):
        p = pair_prob(k, N, Nx, Ny)
        
        p_dist = []
        for k in range(0, nys+1):
            term1 = math.factorial(nys) / (math.factorial(k) * math.factorial(nys - k)) 
            term2 = p**k * (1 - p)**(nys-k)
            p_dist.append(term1 * term2)
        p_dist_ks.append(p_dist)
        
    return p_dist_ks
        
