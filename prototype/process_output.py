# Some functions for processing results - should be put somewhere else
import numpy as np
def get_sig_links(scores,random_scores,p_threshold=0.95,direction = 'greater'):
    # gcfs are columns
    n_spec,n_gcf = scores.shape
    sig_links = np.zeros((n_spec,n_gcf))
    percs = []
    for gpos in range(n_gcf):
        perc = np.percentile(random_scores[:,gpos],p_threshold)
        for i,s in enumerate(scores[:,gpos]):
            if direction == 'greater':
                if s>= perc:
                    sig_links[i,gpos] = 1
            else:
                if s<= perc:
                    sig_links[i,gpos] = 1
    return sig_links


def get_sig_spec(data_link,sig_links,scores,gcf_pos):
    col = sig_links[:,gcf_pos] # get the column
    sig_pos = np.where(col == 1)[0]
    print(sig_pos) 
    orig_ids = []
    for sp in sig_pos:
        orig_ids.append((int(data_link.mapping_spec.iloc[sp]["original spec-id"]),scores[sp,gcf_pos]))
    orig_ids.sort(key=lambda x: x[1],reverse = True)
    return orig_ids