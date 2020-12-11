#!/usr/bin/env python
# coding: utf-8

# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[1]:


from cdk_pywrapper.cdk_pywrapper import Compound


# In[ ]:





# In[ ]:





# In[2]:


from pyteomics import mgf
import csv


# In[3]:


# Convert MIBiG SMILES to InChI keys
mibig_inchi_keys = {}
failing_smiles = []
cnt = 0
with open('compunds_structures_2.0.csv', 'r') as f:
    for l in csv.reader(f):
        if cnt is 0:
            cnt += 1
            continue
        if len(l) is 0:
            continue
        mibig_id, compound_name, smiles, pubchem_id = l
        if smiles == '':
            continue
        # print(smiles)
        mibig_compound = Compound(compound_string=smiles, identifier_type='smiles')
        try:
            mibig_inchi_key = mibig_compound.get_inchi_key()
        except:
            failing_smiles.append(l)
            continue
        mibig_inchi_first_block = mibig_inchi_key.split('-')[0]
        if mibig_inchi_first_block in mibig_inchi_keys:
            mibig_inchi_keys[mibig_inchi_first_block].append(((mibig_id, compound_name, pubchem_id), mibig_inchi_key, smiles))#.split('_')[0])
        else:
            mibig_inchi_keys[mibig_inchi_first_block] = [((mibig_id, compound_name, pubchem_id), mibig_inchi_key, smiles)]#.split('_'[0])]


# In[6]:


# Read InChI keys from spectra
matches = []
gnps_inchi_keys = set()
gnps_count = 0
failing_smiles = []
# for gnps_spectrum in mgf.read('/home/grimur/gnps-mibig/gnps_with_structure.mgf', use_index=False):
for gnps_spectrum in mgf.read('/home/grimur/iokr/data/mibig/gnps-mibig/gnps_with_structure.mgf', use_index=False):
    gnps_inchi_key = gnps_spectrum['params']['inchikey']
    gnps_smiles = gnps_spectrum['params']['smiles']
    try:
        compound = Compound(compound_string=gnps_smiles, identifier_type='smiles')
    except:
        print('Smiles error:')
        print(gnps_smiles)
        failing_smiles.append((gnps_count, gnps_spectrum['params']['filename'], gnps_smiles))
    gnps_count += 1
    print('processing {}'.format(gnps_inchi_key))
    gnps_inchi_first_block = gnps_inchi_key.split('-')[0]
    gnps_inchi_keys.add(gnps_inchi_first_block)
    mibig_id = mibig_inchi_keys.get(gnps_inchi_first_block, None)
    if mibig_id is not None:
        matches.append((gnps_spectrum, mibig_id))


# In[7]:


print('matches: {}'.format(len(matches)))
print('matched mibig ids: {}'.format(len(set([x[1][0][0] for x in matches]))))
print('gnps spectra: {}'.format(gnps_count))
print('gnps keys: {}'.format(len(gnps_inchi_keys)))
print('total mibig ids: {}'.format(len(set().union(*mibig_inchi_keys.values()))))
print('mibig keys: {}'.format(len(mibig_inchi_keys.keys())))


# In[51]:


#print(set([x[1][0][0] for x in matches]))


# In[52]:


# with open('failing_smiles.csv', 'w') as f:
#     for i in failing_smiles:
#         f.write(('%s,%s,%s' % i) + '\n')



# In[6]:


import csv


# In[7]:


# Write a list of the matched MIBiG entries / spectra
with open('matched_mibig_gnps_update.csv', 'w') as f:
    fwriter = csv.writer(f)
    fwriter.writerow(['# mgf_spectrum_id',
                'mgf_inchikey',
                'mgf_smiles',
                'mibig_id',
                'mibig_name',
                'mibig_inchi',
                'mibig_smiles'])
    for spectrum, mibig_entries in matches:
        spec_id = spectrum['params']['spectrumid']
        spec_inchi = spectrum['params']['inchikey']
        spec_smiles = spectrum['params']['smiles']
        for mibig_ids, mibig_inchi, mibig_smiles in mibig_entries:
            mibig_id, mibig_name, ext_id = mibig_ids
            output_list = [spec_id, spec_inchi, spec_smiles, mibig_id, mibig_name, mibig_inchi, mibig_smiles]
            fwriter.writerow(output_list)
            #print(output_string)
            


# In[8]:


# Write a MGF file with the matched spectra
mgf.write([x[0] for x in matches], 'matched_mibig_gnps_update.mgf')


# In[14]:


matches[0][0]['params']


# In[24]:


matches_with_bgc_id = []
for ms, bgc in matches:
    bgc_id = bgc[0][0][0]
    spectrum_id = ms['params']['spectrumid']
    compound_id = '.'.join((bgc_id, spectrum_id))
    ms['params']['BGCID'] = compound_id
    matches_with_bgc_id.append(ms)
    
    


# In[25]:


mgf.write(matches_with_bgc_id, 'matched_mibig_gnps_update_mibigid.mgf')


# In[28]:





# In[ ]:




