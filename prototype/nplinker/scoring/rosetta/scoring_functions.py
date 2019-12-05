from __future__ import print_function

def fast_cosine_shift(spectrum1,spectrum2,tol,min_match):
    if spectrum1.n_peaks == 0 or spectrum2.n_peaks == 0:
        return 0.0,[]

    spec1 = spectrum1.normalised_peaks
    spec2 = spectrum2.normalised_peaks

    zero_pairs = find_pairs(spec1,spec2,tol,shift=0.0)

    shift = spectrum1.parent_mz - spectrum2.parent_mz

    nonzero_pairs = find_pairs(spec1,spec2,tol,shift = shift)

    matching_pairs = zero_pairs + nonzero_pairs

    matching_pairs = sorted(matching_pairs,key = lambda x: x[2], reverse = True)

    used1 = set()
    used2 = set()
    score = 0.0
    used_matches = []
    for m in matching_pairs:
        if not m[0] in used1 and not m[1] in used2:
            score += m[2]
            used1.add(m[0])
            used2.add(m[1])
            used_matches.append(m)
    if len(used_matches) < min_match:
        score = 0.0
    return score,used_matches


def find_pairs(spec1,spec2,tol,shift=0):
    matching_pairs = []
    spec2lowpos = 0
    spec2length = len(spec2)
    
    for idx,(mz,intensity) in enumerate(spec1):
        # do we need to increase the lower idx?
        while spec2lowpos < spec2length and spec2[spec2lowpos][0] + shift < mz - tol:
            spec2lowpos += 1
        if spec2lowpos == spec2length:
            break
        spec2pos = spec2lowpos
        while(spec2pos < spec2length and spec2[spec2pos][0] + shift < mz + tol):
            matching_pairs.append((idx,spec2pos,intensity*spec2[spec2pos][1]))
            spec2pos += 1
        
    return matching_pairs    

def fast_cosine(spectrum1,spectrum2,tol,min_match):
    # spec 1 and spec 2 have to be sorted by mz
    if spectrum1.n_peaks == 0 or spectrum2.n_peaks == 0:
        return 0.0,[]
    # find all the matching pairs
    
    spec1 = spectrum1.normalised_peaks
    spec2 = spectrum2.normalised_peaks
    
    matching_pairs = find_pairs(spec1,spec2,tol,shift = 0.0)
    
        
        
    matching_pairs = sorted(matching_pairs,key = lambda x:x[2],reverse = True)
    used1 = set()
    used2 = set()
    score = 0.0
    used_matches = []
    for m in matching_pairs:
        if not m[0] in used1 and not m[1] in used2:
            score += m[2]
            used1.add(m[0])
            used2.add(m[1])
            used_matches.append(m)
    if len(used_matches) < min_match:
        score = 0.0
    return score,used_matches

def comp_scores(spectra,file_scan,similarity_function,similarity_tolerance,min_match):
    # a method for testing -- just computes scores between a bunch of scans
    specs = []
    for file_name,scan_number in file_scan:
        specs.append(filter(lambda x: x.file_name == file_name and x.scan_number == scan_number,spectra)[0])

    for i in range(len(file_scan)-1):
        for j in range(i+1,len(file_scan)):
            (f,s) = file_scan[i]
            spec = specs[i]
            (f2,s2) = file_scan[j]
            spec2 = specs[j]
            sc,_ = similarity_function(spec,spec2,similarity_tolerance,min_match)
            print("{},{} <-> {},{} = {}".format(f,s,f2,s2,sc))



