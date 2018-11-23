# functions

import numpy as np

def calc_correlation_matrix(M_type1_cond, M_type2_cond):
    """ Calculate correlation matrices from co-occurence matrices
    Input:
    M_type1_cond(x,y) is 1 if type1_x IS observed under condition_y
    M_type1_cond(x,y) is 0 if type1_x IS NOT observed under condition_y
    
    Outputs three correlation matrices:
    M_type1_type2(x,y) --- number of conditions where type1_x and type2_y co-occur
    M_type1_nottype2(x,y) --- number of conditions where type1_x and NOT-type2_y co-occur
    M_nottype1_type2(x,y) --- number of conditions where NOT-type1_x and type2_y co-occur 
    """
       
    dim1 = M_type1_cond.shape[0]
    dim2 = M_type2_cond.shape[0]
    M_type1_type2 = np.zeros((dim1, dim2))
    M_type1_nottype2 = np.zeros((dim1, dim2))
    M_nottype1_type2 = np.zeros((dim1, dim2))
    
    for i in range(0, dim2):
        # show progress:
        if (i+1) % 100 == 0 or i == dim2-1:  
            print('\r', ' Done ', i+1, ' of ', dim2, ' type2s', end="")
            
        # look where type1 AND type2 are present:
        updates = M_type1_cond[:, np.where(M_type2_cond[i,:] > 0)[0]]
        updates = np.sum(updates, axis=1)
        M_type1_type2[:,i] = updates

        # look where type2 i is NOT present:
        updates = M_type1_cond[:, np.where(M_type2_cond[i,:] == 0)[0]]
        updates = np.sum(updates, axis=1)
        M_type1_nottype2[:,i] = updates
    print("")
    
    for i in range(0, dim1):
        # show progress:
        if (i+1) % 100 == 0 or i == dim1-1:
            print('\r', ' Done ', i+1, ' of ', dim1, ' type1s.', end="")

        # look where type1 i is NOT present:
        updates = M_type2_cond[:, np.where(M_type1_cond[i,:] == 0)[0]]
        updates = np.sum(updates, axis=1)
        M_nottype1_type2[i,:] = updates
    print("")
    
    return M_type1_type2, M_type1_nottype2, M_nottype1_type2 



def calc_likelihood_matrix(M_type1_cond, M_type2_cond, 
                           M_type1_type2, M_type1_nottype2, M_nottype1_type2):
    """ Calculate correlation matrices from co-occurence matrices
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