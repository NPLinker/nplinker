# functions

import numpy as np
import math

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
        

def pair_prob_approx(P_str, XG, Ny, hits):
    """
    Calculate probability of finding 'k' hits between Gx and Sy.
    
    Parameters
    ----------
    P_str: numpy array
        Probabilities for finding a spectrum in the a certain strain.
        Usually this can be set to num-of-spectra-in-strain / num-of-spectra-in-all-strains
    XG: list
        List of ids of strains where the GCF of interest occurs.
    Ny: int
        Number of strains that contain the spectrum of interest.
    hits: int
        number of hits
    """ 
    
    Nx = len(XG)
    Nstr = len(P_str)
    
    # Check Nx, Ny, hits
    if (hits > Nx) or (hits > Ny):
        print("Given number of 'hits' must be <= Nx and <= Ny.")

    p_hit_mean = np.sum(P_str[XG])/Nx
    p_nohit_mean = (1 - np.sum(P_str[XG]))/(Nstr - Nx)
    p_mean = (hits * p_hit_mean + (Ny-hits) * p_nohit_mean)/Ny
    
    # Calculate product of all hits
    prod0 = 1
    for i in range(hits):
        prod0 = prod0 * p_hit_mean * (Nx -i)

    # Calculate product of all non-hits
    prod1 = 1
    for i in range(Ny - hits):
        prod1 = prod1 * ((Nstr- Nx - i)/Nstr)

    # Calculate product of probability updates 
    # (fewer accessible elements lead to increasing probabilities)
    prod2 = 1
    for j in range(Ny):
        prod2 = prod2 * (1/(1 - j*p_mean))  
    
    return np.sum(math.factorial(Ny)/(math.factorial(hits) * math.factorial(Ny - hits)) * prod0 * prod1 * prod2)



def link_prob(P_str, XGS, Nx, Ny, Nstr):
    """
    Calculate probability of finding a set of *specific* hits between Gx and Sy.
    This means: the probability to find hits only in all the strains that form the set XGS
    
    Parameters
    ----------
    P_str: numpy array
        Probabilities for finding a spectrum in the a certain strain.
        Usually this can be set to num-of-spectra-in-strain / num-of-spectra-in-all-strains
    XGS: list
        List of ids of strains where GCF and a spectrum of interest co-occurs.    
    Nx: int
        Number of strains that contain the GCF of interest.
    Ny: int
        Number of strains that contain the spectrum of interest.
    Nstr: int
        Number of strains.
    """ 
    p_mean = 1/Nstr
    hits = len(XGS)
    
    # Calculate product of all non-hits
    prod1 = 1
    for i in range(Ny - hits):
        prod1 = prod1 * ((Nstr- Nx - i)/Nstr)

    # Calculate product of probability updates 
    # (fewer accessible elements lead to increasing probabilities)
    prod2 = 1
    for j in range(Ny):
        prod2 = prod2 * (1/(1 - j*p_mean))  
        
    return math.factorial(Ny)/math.factorial(Ny-hits) * np.prod(P_str[XGS]) * prod1 * prod2



def pair_prob_hg(k, N, Nx, Ny):
    """
    Calculate the probability to draw k times type(Ny) out of N elements (whereof Ny type(Ny)s),
    when drawing Nx times in total.
    Same as hypergemoetric distribution
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
        p = pair_prob_hg(k, N, Nx, Ny)
        
        p_dist = []
        for k in range(0, nys+1):
            term1 = math.factorial(nys) / (math.factorial(k) * math.factorial(nys - k)) 
            term2 = p**k * (1 - p)**(nys-k)
            p_dist.append(term1 * term2)
        p_dist_ks.append(p_dist)
        
    return p_dist_ks
        

def pair_prob(P_str, XG, Ny, hits):
    """
    Calculate probability of finding 'k' hits between Gx and Sy.
         
    CAREFUL: for larger Nx, Ny, Nstr this quickly becomes *VERY* slow (many, many permutations...)
    --> better use pair_prob_approx instead
    
    Parameters
    ----------
    P_str: numpy array
        Probabilities for finding a spectrum in the a certain strain.
        Usually this can be set to num-of-spectra-in-strain / num-of-spectra-in-all-strains
    XG: list
        List of ids of strains where the GCF of interest occurs.
    Ny: int
        Number of strains that contain the spectrum of interest.
    hits: int
        number of hits
    """ 
    
    Nx = len(XG)
    Nstr = len(P_str)
    
    # Check Nx, Ny, hits
    if (hits > Nx) or (hits > Ny):
        print("Given number of 'hits' must be <= Nx and <= Ny.")

    # Calculate all unique permutations:        
    state0 = [1]*hits + [0]*(Nx-hits)
    states = np.array(list(permutation_unique(state0)))
    
    # Calculate the product of all probabilties accross all permutations
    P_states = states*P_str[XG]
    prods = np.prod(P_states + np.abs(states-1), axis=1)
    del P_states
    del states
    p_mean = 1/Nstr
    
    # Calculate product of all non-hits
    prod1 = 1
    for i in range(Ny - hits):
        prod1 = prod1 * ((Nstr- Nx - i)/Nstr)

    # Calculate product of probability updates 
    # (fewer accessible elements lead to increasing probabilities)
    prod2 = 1
    for j in range(Ny):
        prod2 = prod2 * (1/(1 - j*p_mean))  
    
    return np.sum(math.factorial(Ny)/math.factorial(Ny-hits) * prods * prod1 * prod2)


# method to calculate unique permutations:
class unique_element:
    def __init__(self,value,occurrences):
        self.value = value
        self.occurrences = occurrences

def permutation_unique(elements):
    """
    Derive unique permutations of elements (list)
    """
    eset = set(elements)
    listunique = [unique_element(i, elements.count(i)) for i in eset]
    num_elements = len(elements)
    return permutation_unique_helper(listunique, [0]*num_elements, num_elements-1)

def permutation_unique_helper(listunique, result_list, d):
    """
    Helper function to derive unique permutations of elements (list)
    """
    if d < 0:
        yield tuple(result_list)
    else:
        for i in listunique:
            if i.occurrences > 0:
                result_list[d]=i.value
                i.occurrences-=1
                for g in  permutation_unique_helper(listunique, result_list, d-1):
                    yield g
                i.occurrences+=1
