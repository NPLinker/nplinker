# Methods to find correlations between spectra/molecular families and 
# gene clusters/families (BGCs/GCFs)
#
# (still at very much protoype/exploration stage)
#
# Naming:
# M_*   stands for a matrix format
# map_* stands for a simple mapping lookup table
# spec  stands for spectrum
# fam   stands for molecular family

# import packages
import numpy as np
from collections import Counter
import pandas as pd

# import packages for plotting
from matplotlib import pyplot as plt
import seaborn as sns


class DataLinks(object):
    def __init__(self):
        # matrices that store co-occurences with respect to strains
        self.M_strain_gcf = []
        self.M_strain_spec = []
        self.M_strain_fam = []

        # mappings (lookup lists to map between different ids and categories
        self.map_spec_id = []
        self.map_spec_fam = []
        self.map_gcf_class = []
        self.map_strain_name = []

        # correlation matrices for spectra <-> GCFs
        self.M_spec_gcf = []  # = 1 where spec_x and gcf_y co_occure
        self.M_spec_notgcf = []  # = 1 where spec_x and NO gcf_y co_occure
        self.M_notspec_gcf = []  # = 1 where NO spec_x and gcf_y co_occure
        # and the same for mol.families <-> GCFs
        self.M_fam_gcf = []
        self.M_fam_notgcf = []
        self.M_notfam_gcf = []
        self.strain_fam_labels = []  # labels for strain-family matrix


    def load_data(self, spectra, gcf_list, strain_list):
        # load data from spectra, GCFs, and strains
        print("Create mappings between spectra, gcfs, and strains.")
        self.mappings_spec(spectra)
        self.mappings_gcf(gcf_list)
        print("Create co-occurence matrices: spectra<->strains + and gcfs<->strains.")
        self.matrix_strain_gcf(gcf_list, strain_list)
        self.matrix_strain_spec(spectra, strain_list)

      
    def find_correlations(self, include_singletons=False):
        # collect correlations/ co-occurences
        print("Create correlation matrices: spectra<->gcfs.")
        self.correlation_matrices(type='spec-gcf')
        print("Create correlation matrices: mol-families<->gcfs.")
        self.data_family_mapping(include_singletons=include_singletons)
        self.correlation_matrices(type='fam-gcf')


    def mappings_spec(self, spectra):
        # Collect most import mapping tables from input data
        self.map_spec_id = np.zeros((len(spectra),1))
        self.map_spec_fam = np.zeros((len(spectra),1))

        if isinstance(spectra[0].family, str):  
            for i, spectrum in enumerate(spectra):
                self.map_spec_fam[i,0] = spectrum.family
                self.map_spec_id[i,0] = spectrum.spectrum_id
        elif isinstance(spectra[0].family.family_id, str): # assume make_families was run
            for i, spectrum in enumerate(spectra):
                self.map_spec_fam[i,0] = spectrum.family.family_id
                self.map_spec_id[i,0] = spectrum.spectrum_id
        else:
            print("No proper family-id found in spectra.")


    def mappings_gcf(self, gcf_list):
        # find classes of gene cluster (nrps, pksi etc.)
        # collect most likely class (most occuring name, preferentially not "Others")
        # additional score shows fraction of chosen class among all given ones
        bigscape_bestguess = []
        for i, gcf in enumerate(gcf_list):
            bigscape_class = []
            for m in range(0, len(gcf_list[i].bgc_list)):
                bigscape_class.append(gcf_list[i].bgc_list[m].bigscape_class)
                class_counter = Counter(bigscape_class)
            
            # try not to select "Others"
            if class_counter.most_common(1)[0][0] == "Others" and class_counter.most_common(1)[0][1] < len(bigscape_class): 
                bigscape_bestguess.append([class_counter.most_common(2)[1][0], class_counter.most_common(2)[1][1]/len(bigscape_class)])
            else:
                bigscape_bestguess.append([class_counter.most_common(1)[0][0], class_counter.most_common(1)[0][1]/len(bigscape_class)])

        self.map_gcf_class = bigscape_bestguess


    def matrix_strain_gcf(self, gcf_list, strain_list):
        # Collect co-ocurences in M_strain_spec matrix
        M_strain_gcf = np.zeros((len(gcf_list), len(strain_list)))

        for i, strain  in enumerate(strain_list):
            for m, gcf in enumerate(gcf_list):
                in_gcf = gcf.has_strain(strain)
                if in_gcf:
                    M_strain_gcf[m,i] += 1

        self.M_strain_gcf = M_strain_gcf


    def matrix_strain_spec(self, spectra, strain_list):
        # Collect co-ocurences in M_strains_specs matrix
        M_strain_spec = np.zeros((len(spectra), len(strain_list)))

        # find strain data in spectra metadata
        metadata_category, metadata_value = zip(*list(spectra[0].metadata.items()))
        strain_list = [str(x) for x in strain_list]
        strain_meta_num = []
        for num, category in enumerate(metadata_category):
            if category in strain_list:
                strain_meta_num.append(num)
        # TODO add test to see if all strain names are found in metadata

        for i,spectrum in enumerate(spectra):
            metadata_category, metadata_value = zip(*list(spectra[i].metadata.items()))
            M_strain_spec[i, 0:len(strain_meta_num)] = [metadata_value[x] for x in strain_meta_num]

        strain_spec_labels = [metadata_category[x] for x in strain_meta_num]  # strain names

        self.M_strain_spec = M_strain_spec
        self.map_strain_name = strain_spec_labels


    def data_family_mapping(self, include_singletons=False):
        # Create M_strain_fam matrix that gives co-occurences between mol. families and strains
        # matrix dimensions are: number of families  x  number of strains

        family_ids = np.unique(self.map_spec_fam) # get unique family ids

        if include_singletons and np.where(self.map_spec_fam == -1)[0].shape[0] > 0: 
            num_of_singletons = np.where(self.map_spec_fam == -1)[0].shape[0]
            num_unique_fams = num_of_singletons + len(family_ids) - 1
        else:
            num_of_singletons = 0
            num_unique_fams = len(family_ids)
   
        M_strain_fam = np.zeros((num_unique_fams, self.M_strain_spec.shape[1]))
        strain_fam_labels = []

        if num_of_singletons > 0:  # if singletons exist + included
            M_strain_fam[(num_unique_fams-num_of_singletons):,
                         :] = self.M_strain_spec[np.where(self.map_spec_fam[:,0] == -1)[0],:]     
            
        # go through families (except singletons) and collect member strain occurences    
        for i, fam_id in enumerate(family_ids[np.where(family_ids != -1)].astype(int)):
            M_strain_fam[i,:] = np.sum(self.M_strain_spec[np.where(self.map_spec_fam[:,0] == fam_id),:], axis=1)
            strain_fam_labels.append(fam_id)
            
        # only looking for co-occurence, hence only 1 or 0
        M_strain_fam[M_strain_fam>1] = 1
        strain_fam_labels.append([-1] * num_of_singletons)
        
        self.M_strain_fam = M_strain_fam
        self.strain_fam_labels = strain_fam_labels


    def correlation_matrices(self, type='spec-gcf'):
        # collect co-ocurences:
        # IF type='spec-gcf':  spectra and GCFS in M_spec_gcf matrix
        # IF type='fam-gcf':  mol.families and GCFS in M_fam_gcf matrix

        if type == 'spec-gcf':
            M_strain_type1 = self.M_strain_spec
        elif type == 'fam-gcf':
            M_strain_type1 = self.M_strain_fam
        elif type == 'spec-bgc' or type == 'fam-bgc':
            print("Given types are not yet supported... ")
        else:
            print("Wrong correlation 'type' given.")
            print("Must be one of 'spec-gcf', 'fam-gcf', ...")

        print("Calculating correlation matrices of type: ", type)
        dim1 = M_strain_type1.shape[0]
        dim2 = self.M_strain_gcf.shape[0]

        M_type1_gcf = np.zeros((dim1, dim2))
        M_type1_notgcf = np.zeros((dim1, dim2))
        M_nottype1_gcf = np.zeros((dim1, dim2))

        for i in range(0, dim2):
            if (i+1) % 100 == 0 or i == dim2-1:  # show progress
                print('\r', ' Done ', i+1, ' of ', dim2, ' GCFs.', end="")
                
            updates = M_strain_type1[:, np.where(self.M_strain_gcf[i,:] > 0)[0]]
            updates = np.sum(updates, axis=1)

            M_type1_gcf[:,i] = updates
            # look where GCF i is NOT present:
            updates = M_strain_type1[:, np.where(self.M_strain_gcf[i,:] == 0)[0]]
            updates = np.sum(updates, axis=1)
            M_type1_notgcf[:,i] = updates
        print("")

        for i in range(0, dim1):
            if (i+1) % 100 == 0 or i == dim1-1:  # show progress
                print('\r', ' Done ', i+1, ' of ', dim1, ' spectra/families.', end="")
            # look where spectrum or family i is NOT present:
            updates = self.M_strain_gcf[:, np.where(M_strain_type1[i,:] == 0)[0]]
            updates = np.sum(updates, axis=1)
            M_nottype1_gcf[i,:] = updates

        # return results:
        if type == 'spec-gcf':
            self.M_spec_gcf = M_type1_gcf
            self.M_spec_notgcf = M_type1_notgcf
            self.M_notspec_gcf = M_nottype1_gcf
        elif type == 'fam-gcf':
            self.M_fam_gcf = M_type1_gcf
            self.M_fam_notgcf = M_type1_notgcf
            self.M_notfam_gcf = M_nottype1_gcf
        else:
            print("No correct correlation matrix was created.")
        print("")

    # class data_links OUTPUT functions
    # TODO add output functions (e.g. to search for mappings of individual specs, gcfs etc.)



class LinkProbability(object):
    def __init__(self):
        # matrices that store probabilities of co-occurence
        # P_spec_givengcf contains probabilities P(spec_x|gcf_y),
        # which is the probability of finding spec_x given there is a gcf_y

        # Probabilities reflecting co-occurences spectra <-> GCFs
        self.P_spec_given_gcf = []
        self.P_spec_not_gcf = []
        self.P_gcf_given_spec = []
        self.P_gcf_not_spec = []
        # Probabilities reflecting co-occurences mol.families <-> GCFs
        self.P_fam_given_gcf = []
        self.P_fam_not_gcf = []
        self.P_gcf_given_fam = []
        self.P_gcf_not_fam = []
        # selection of candidates
        self.link_candidates_gcf_spec = []
        self.link_candidates_gcf_fam = []
#        self.candidates_fam_gcf = []
#        self.gcf_fam_competition = []
#        self.fam_gcf_winners = []
#        self.fam_gcf_winnersscores = []


    def correlation_probabilities(self, data_links, type='spec-gcf'):
        # Calulate probabilites that reflect co-occurences in data
        # IF type='spec-gcf':  
        # P(GCF_x | spec_y), P(spec_y | GCF_x), 
        # P(GCF_x | not spec_y), P(spec_y | not GCF_x) 
        # IF type='fam-gcf':  
        # P(GCF_x | fam_y), P(fam_y | GCF_x), 
        # P(GCF_x | not fam_y), P(fam_y | not GCF_x) 

        if type == 'spec-gcf':
            M_type1_gcf = data_links.M_spec_gcf 
            M_type1_notgcf = data_links.M_spec_notgcf
            M_nottype1_gcf = data_links.M_notspec_gcf
            M_strain_type1 = data_links.M_strain_spec
        elif type == 'fam-gcf':
            M_type1_gcf = data_links.M_fam_gcf
            M_type1_notgcf = data_links.M_fam_notgcf
            M_nottype1_gcf = data_links.M_notfam_gcf
            M_strain_type1 = data_links.M_strain_fam
        elif type == 'spec-bgc' or type == 'fam-bgc':
            print("Given types are not yet supported... ")
        else:
            print("Wrong correlation 'type' given.")
            print("Must be one of 'spec-gcf', 'fam-gcf'...")

        print("Calculating probability matrices of type: ", type)
        dim1, dim2 = M_type1_gcf.shape
        num_strains = data_links.M_strain_gcf.shape[1]

        P_gcf_given_type1 = np.zeros((dim1, dim2))
        P_gcf_not_type1 = np.zeros((dim1, dim2))
        P_type1_given_gcf = np.zeros((dim1, dim2))
        P_type1_not_gcf = np.zeros((dim1, dim2))

        # Calculate P_gcf_given_type1 matrix
        P_sum_spec = np.sum(M_strain_type1, axis=1)
        P_sum_spec[P_sum_spec < 1] = 1 #avoid later division by 0
        P_gcf_given_type1 = M_type1_gcf/np.tile(P_sum_spec, 
                                                (dim2, 1)).transpose(1, 0)

        # Calculate P_type1_given_gcf matrix
        P_sum_gcf = np.sum(data_links.M_strain_gcf, axis=1)
        P_sum_gcf[P_sum_gcf < 1] = 1 #avoid later division by 0
        P_type1_given_gcf = M_type1_gcf/np.tile(P_sum_gcf, (dim1, 1))

        # Calculate P_gcf_not_type1 matrix
        P_sum_not_spec = num_strains - P_sum_spec
        P_gcf_not_type1 = M_nottype1_gcf/np.tile(P_sum_not_spec,
                                                 (dim2, 1)).transpose(1, 0)

        # Calculate P_type1_not_gcf matrix
        P_sum_not_gcf = num_strains - P_sum_gcf  
        P_type1_not_gcf = M_type1_notgcf/np.tile(P_sum_not_gcf, (dim1, 1))

        if type == 'spec-gcf':
            self.P_gcf_given_spec = P_gcf_given_type1
            self.P_gcf_not_spec = P_gcf_not_type1
            self.P_spec_given_gcf = P_type1_given_gcf
            self.P_spec_not_gcf = P_type1_not_gcf
        elif type == 'fam-gcf':
            self.P_gcf_given_fam = P_gcf_given_type1
            self.P_gcf_not_fam = P_gcf_not_type1
            self.P_fam_given_gcf = P_type1_given_gcf
            self.P_fam_not_gcf = P_type1_not_gcf
        else:
            print("No correct probability matrices were created.")


    def select_link_candidates(self, data_links, P_cutoff=0.8, score_cutoff=0, type='fam-gcf'):
        # look for potential best candidate for links between
        # IF type='spec-gcf': GCFs and spectra
        # IF type='fam-gcf': GCFs and mol.families
        #
        # Only consider candidates if P_gcf_given_type1 >= P_cutoff

        if type == 'spec-gcf':
            P_gcf_given_type1 = self.P_gcf_given_spec
            P_gcf_not_type1 = self.P_gcf_not_spec
            P_type1_given_gcf = self.P_spec_given_gcf
            P_type1_not_gcf = self.P_spec_not_gcf
            M_type1_gcf = data_links.M_spec_gcf
            index_names = ["spectrum_id", "GCF id", "P(gcf|spec)", "P(spec|gcf)", 
                             "P(gcf|not spec)", "P(spec|not gcf)", 
                             "co-occur in # strains", "link score"]
        elif type == 'fam-gcf':
            P_gcf_given_type1 = self.P_gcf_given_fam
            P_gcf_not_type1 = self.P_gcf_not_fam
            P_type1_given_gcf = self.P_fam_given_gcf
            P_type1_not_gcf = self.P_fam_not_gcf
            M_type1_gcf = data_links.M_fam_gcf
            index_names = ["family_id", "GCF id", "P(gcf|fam)", "P(fam|gcf)", 
                             "P(gcf|not fam)", "P(fam|not gcf)", 
                             "co-occur in # strains", "link score"]
        elif type == 'spec-bgc' or type == 'fam-bgc':
            print("Given types are not yet supported... ")
        else:
            print("Wrong correlation 'type' given.")
            print("Must be one of 'spec-gcf', 'fam-gcf'...")

        dim1, dim2 = P_gcf_given_type1.shape

        # ROUND 1: select candidates with P_gcf_given_spec >= P_cutoff:
        candidate_ids = np.where(P_gcf_given_type1[:,:] >= P_cutoff)
        link_candidates = np.zeros((8, candidate_ids[0].shape[0]))
        link_candidates[0,:] = candidate_ids[0].astype(int)  # spectrum/fam number
        link_candidates[1,:] = candidate_ids[1].astype(int)  # gcf id
        link_candidates[2,:] = P_gcf_given_type1[candidate_ids]
        link_candidates[3,:] = P_type1_given_gcf[candidate_ids]
        link_candidates[4,:] = P_gcf_not_type1[candidate_ids]
        link_candidates[5,:] = P_type1_not_gcf[candidate_ids]
        link_candidates[6,:] = M_type1_gcf[candidate_ids].astype(int)

        # The next is the actual SCORING. Just a starting point...
        # now it is:
        # P_gcf_given_type1 * (1-P_type1_not_gcf) weighted by no. of strains they co-occure
        link_candidates[7,:] = link_candidates[2,:] * (1 - link_candidates[5,:]) \
                            *(1-np.exp(-0.5 * link_candidates[6,:]))
        print(link_candidates.shape[1], " candidates selected with ", 
              index_names[2], " >= ", P_cutoff, " .")

        # ROUND 2: select based on score >= score_cutoff
        link_candidates = link_candidates[:, link_candidates[7,:] >= score_cutoff]
        print(link_candidates.shape[1], " candidates selected with ", 
              index_names[2], " >= ", P_cutoff, " and a link score >= ", score_cutoff, ".")
        
        # Transform into pandas Dataframe (to avoid index confusions):
        link_candidates_pd = pd.DataFrame(link_candidates.transpose(1,0), columns = index_names)
        
        # add other potentially relevant knowdledge
        # If this will grow to more collected information -> create separate function
        bgc_class = [] 
        for i in link_candidates_pd["GCF id"].astype(int):
            bgc_class.append(data_links.map_gcf_class[i][0])
        link_candidates_pd["BGC class"] = bgc_class

        # return results
        if type == 'spec-gcf':
            self.link_candidates_gcf_spec = link_candidates_pd
        elif type == 'fam-gcf':   
            self.link_candidates_gcf_fam = link_candidates_pd
        else:
            print("No candidate selection was created.")

    
    def plot_candidates(self, P_cutoff=0.8, score_cutoff=0, type='fam-gcf'):
        # plot best rated correlations between gcfs and spectra/families
        # plot in form of seaborn clustermap
        
        if type == 'spec-gcf':
            link_candidates = self.link_candidates_gcf_spec
            selected_ids = np.where((link_candidates["P(gcf|spec)"] > P_cutoff) & 
                        (link_candidates["link score"] > score_cutoff))[0]
        elif type == 'fam-gcf':   
            link_candidates = self.link_candidates_gcf_fam
            selected_ids = np.where((link_candidates["P(gcf|fam)"] > P_cutoff) & 
                        (link_candidates["link score"] > score_cutoff))[0]
        else:
            print("Wrong correlation 'type' given.")   
    
        mapping_fams = np.unique(link_candidates.iloc[selected_ids, 0])
        mapping_gcfs = np.unique(link_candidates.iloc[selected_ids, 1])
        unique_fams = len(mapping_fams)
        unique_gcfs = len(mapping_gcfs)
    
        M_links = np.zeros((unique_fams, unique_gcfs))
    
        # define colors for different BGC classes
        bigscape_classes_dict = {"Others": "C0", 
                                 "NRPS": "C1", 
                                 "PKS-NRP_Hybrids": "C2",
                                 "PKSother": "C3", 
                                 "PKSI": "C4", 
                                 "RiPPs": "C5", 
                                 "Saccharides": "C6",
                                 "Terpene": "C7"}
        col_colors = []
    
        # Create matrix with relevant link scores
        # TODO replace by better numpy method...
        for i in range(0, link_candidates.shape[0]):
            x = np.where(mapping_fams == link_candidates.iloc[i, 0])[0]
            y = np.where(mapping_gcfs == link_candidates.iloc[i, 1])[0]
            M_links[x, y] = link_candidates["link score"][i] 
    
        # make pandas dataframe from numpy array
        M_links = pd.DataFrame(M_links, 
                               index = mapping_fams.astype(int), 
                               columns = mapping_gcfs.astype(int))
        if type == 'spec-gcf': 
            M_links.index.name = 'spectrum number'
        elif type == 'fam-gcf':
            M_links.index.name = 'molecular family number'
        M_links.columns.name = 'gene cluster family (GCF)'
    
        # add color label representing gene cluster class
        for bgc_class in link_candidates["BGC class"]:
            if bgc_class is None:
                col_colors.append((0,0,0))
            else:
                col_colors.append(bigscape_classes_dict[bgc_class])
    
    #    bgc_type_colors = pd.Series(mapping_gcfs.astype(int), index=M_links.columns).map(col_colors)
        graph = sns.clustermap(M_links, metric="correlation", method="weighted", 
                               cmap="Reds", #sns.cubehelix_palette(8, start=0, rot=.3, dark=0, light=1),
                               vmin=0, vmax=np.max(link_candidates["link score"]), 
                               col_cluster=True, 
                               col_colors=col_colors, 
                               robust=True)
        graph.fig.suptitle('Correlation map') 
    
        # Rotate labels
        plt.setp(graph.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
        plt.setp(graph.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    
        # Make labels smaller
        plt.setp(graph.ax_heatmap.xaxis.get_majorticklabels(), fontsize=7)
        plt.setp(graph.ax_heatmap.yaxis.get_majorticklabels(), fontsize=7)
        
        plt.ylabel("scoring index")
    
        return M_links