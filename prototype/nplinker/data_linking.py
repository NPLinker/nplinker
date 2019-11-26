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
from scipy.stats import hypergeom

# import packages for plotting
# TODO move plotting to separate module?
try:
    from matplotlib import pyplot as plt
    import seaborn as sns
except ImportError:
    print('Warning: plotting functionality will not be available (missing matplotlib and/or seaborn)')

from .genomics import GCF
from .metabolomics import Spectrum
from .metabolomics import MolecularFamily

from .data_linking_functions import calc_correlation_matrix, calc_likelihood_matrix
from .data_linking_functions import pair_prob, pair_prob_hg, link_prob, pair_prob_approx

SCORING_METHODS = ['metcalf', 'likescore', 'hg']

from .logconfig import LogConfig
logger = LogConfig.getLogger(__file__)

class DataLinks(object):
    """ 
    DataLinks collects and structures co-occurence data
    1) Co-occurences of spectra, families, and GCFs with respect to strains
    2) Mappings: Lookup-tables that link different ids and categories
    3) Correlation matrices that show how often spectra/families and GCFs co-occur
    """

    def __init__(self):
        # matrices that store co-occurences with respect to strains
        # values = 1 where gcf/spec/fam occur in strain
        # values = 0 where gcf/spec/fam do not occur in strain
        self.M_gcf_strain = []
        self.M_spec_strain = []
        self.M_fam_strain = []

        # mappings (lookup lists to map between different ids and categories
        self.mapping_spec = pd.DataFrame()
        self.mapping_gcf = pd.DataFrame()
        self.mapping_fam = pd.DataFrame()  # labels for strain-family matrix
        self.mapping_strain = pd.DataFrame()
        self.family_members = []

        # correlation matrices for spectra <-> GCFs
        self.M_spec_gcf = []  # = int: Number of strains where spec_x and gcf_y co_occure
        self.M_spec_notgcf = []  # = int: Number of strains where spec_x and NOT-gcf_y co_occure
        self.M_notspec_gcf = []  # = int: Number of strains where NOT-spec_x and gcf_y co_occure
        # and the same for mol.families <-> GCFs
        self.M_fam_gcf = []
        self.M_fam_notgcf = []
        self.M_notfam_gcf = []


    def get_spec_pos(self,spec_id):
        # get the position in the arrays of a spectrum
        row = self.mapping_spec.loc[self.mapping_spec['original spec-id'] == float(spec_id)]
        return int(row.iloc[0]['spec-id'])

    def get_gcf_pos(self,gcf_id):
        # TODO: fix this so the original ID is present in case of re-ordering
        pass

    def load_data(self, spectra, gcf_list, strain_list):
        # load data from spectra, GCFs, and strains
        logger.debug("Create mappings between spectra, gcfs, and strains.")
        self.collect_mappings_spec(spectra)
        self.collect_mappings_gcf(gcf_list)
        logger.debug("Create co-occurence matrices: spectra<->strains + and gcfs<->strains.")
        self.matrix_strain_gcf(gcf_list, strain_list)
        self.matrix_strain_spec(spectra, strain_list)

    def find_correlations(self, include_singletons=False):
        # collect correlations/ co-occurences
        logger.debug("Create correlation matrices: spectra<->gcfs.")
        self.correlation_matrices(type='spec-gcf')
        logger.debug("Create correlation matrices: mol-families<->gcfs.")
        self.data_family_mapping(include_singletons=include_singletons)
        self.correlation_matrices(type='fam-gcf')


    def collect_mappings_spec(self, spectra):
        # Collect most import mapping tables from input data
        mapping_spec = np.zeros((len(spectra),3))
        mapping_spec[:,0] = np.arange(0,len(spectra))

        if isinstance(spectra[0].family, int):
            for i, spectrum in enumerate(spectra):
                mapping_spec[i,1] = spectrum.id
                mapping_spec[i,2] = spectrum.family
        elif isinstance(spectra[0].family.family_id, int): # assume make_families was run
            for i, spectrum in enumerate(spectra):
                mapping_spec[i,1] = spectrum.id
                mapping_spec[i,2] = spectrum.family.family_id
        else:
            raise Exception("No proper family-id found in spectra.")
            
        # extend mapping tables:
        self.mapping_spec["spec-id"] = mapping_spec[:,0]
        self.mapping_spec["original spec-id"] = mapping_spec[:,1]
        self.mapping_spec["fam-id"] = mapping_spec[:,2]
    
    def collect_mappings_gcf(self, gcf_list):
        """
        Find classes of gene cluster (nrps, pksi etc.)
        collect most likely class (most occuring name, preferentially not "Others")
        additional score shows fraction of chosen class among all given ones
        """
        
        # TODO: not only collect bigclass types but also product preditions
        
        bigscape_bestguess = []
        for i, gcf in enumerate(gcf_list):
            bigscape_class = []
            # for m in range(0, len(gcf_list[i].bgc_list)):
            for i, bgc in enumerate(gcf_list[i].bgcs):
                # bigscape_class.append(gcf_list[i].bgc_list[m].bigscape_class)
                bigscape_class.append(bgc.bigscape_class)
                class_counter = Counter(bigscape_class)
                
            # try not to select "Others":   
            if class_counter.most_common(1)[0][0] == None:
                bigscape_bestguess.append(["Others", 0])
            elif class_counter.most_common(1)[0][0] == "Others" and class_counter.most_common(1)[0][1] < len(bigscape_class): 
                if class_counter.most_common(2)[1][0] == None:
                    bigscape_bestguess.append([class_counter.most_common(1)[0][0], class_counter.most_common(1)[0][1]/len(bigscape_class)])
                else:
                    bigscape_bestguess.append([class_counter.most_common(2)[1][0], class_counter.most_common(2)[1][1]/len(bigscape_class)])
            else:
                bigscape_bestguess.append([class_counter.most_common(1)[0][0], class_counter.most_common(1)[0][1]/len(bigscape_class)])
                
#            if bigscape_bestguess[-1] == None :
#                bigscape_bestguess[-1] = ["Others", 0]

        # extend mapping tables:
        self.mapping_gcf["gcf-id"] = np.arange(0, len(bigscape_bestguess))
        bigscape_guess, bigscape_guessscore = zip(*bigscape_bestguess)
        self.mapping_gcf["bgc class"] = bigscape_guess
        self.mapping_gcf["bgc class score"] = bigscape_guessscore


    def matrix_strain_gcf(self, gcf_list, strain_list):
        # Collect co-ocurences in M_spec_strain matrix
        M_gcf_strain = np.zeros((len(gcf_list), len(strain_list)))

        for i, strain  in enumerate(strain_list):
            for m, gcf in enumerate(gcf_list):
                in_gcf = gcf.has_strain(strain)
                if in_gcf:
                    M_gcf_strain[m,i] = 1

        self.M_gcf_strain = M_gcf_strain
        # extend mapping tables:
        self.mapping_gcf["no of strains"] =  np.sum(self.M_gcf_strain, axis=1)
        self.mapping_strain["no of gcfs"] =  np.sum(self.M_gcf_strain, axis=0)

    def matrix_strain_spec(self, spectra, strain_list):
        # Collect co-ocurences in M_strains_specs matrix

        M_spec_strain = np.zeros((len(spectra), len(strain_list)))
        for i, spectrum in enumerate(spectra):
            for j, s in enumerate(strain_list):
                if spectrum.has_strain(s):
                    M_spec_strain[i,j] = 1
        self.M_spec_strain = M_spec_strain
        
        # extend mapping tables:
        self.mapping_spec["no of strains"] =  np.sum(self.M_spec_strain, axis=1)
        self.mapping_strain["no of spectra"] =  np.sum(self.M_spec_strain, axis=0)
        self.mapping_strain["strain name"] = [str(s) for s in strain_list]

    def data_family_mapping(self, include_singletons=False):
        # Create M_fam_strain matrix that gives co-occurences between mol. families and strains
        # matrix dimensions are: number of families  x  number of strains

        family_ids = np.unique(self.mapping_spec["fam-id"]) # get unique family ids

        if include_singletons and np.where(self.mapping_spec["fam-id"] == -1)[0].shape[0] > 0:
            num_of_singletons = np.where(self.mapping_spec["fam-id"] == -1)[0].shape[0]
            num_unique_fams = num_of_singletons + len(family_ids) - 1
        else:
            num_of_singletons = 0
            num_unique_fams = len(family_ids)

        M_fam_strain = np.zeros((num_unique_fams, self.M_spec_strain.shape[1]))
        strain_fam_labels = []
        strain_fam_index = []
        
        if num_of_singletons > 0:  # if singletons exist + included
            M_fam_strain[(num_unique_fams-num_of_singletons):,
                         :] = self.M_spec_strain[np.where(self.mapping_spec["fam-id"][:,0] == -1)[0],:]

        # go through families (except singletons) and collect member strain occurences  
        self.family_members = []
        for i, fam_id in enumerate(family_ids[np.where(family_ids != -1)].astype(int)):
            family_members = np.where(np.array(self.mapping_spec["fam-id"]) == fam_id)
            self.family_members.append(family_members)
            M_fam_strain[i,:] = np.sum(self.M_spec_strain[family_members,:], axis=1)
            strain_fam_labels.append(fam_id)
            strain_fam_index.append(i)

        add_singleton_entries = -1 in family_ids
        # TODO: i think this breaks stuff below due to mismatches in the number of rows
        # in the dataframes and matrices if there are no -1 family ids. 
        # discovered when trying to write some code to test scoring. is this ever
        # likely to happen with a real dataset??
        if add_singleton_entries:
            strain_fam_labels.append([-1] * num_of_singletons)
            strain_fam_index.append(i+1)
        
        # only looking for co-occurence, hence only 1 or 0
        M_fam_strain[M_fam_strain>1] = 1

        self.M_fam_strain = M_fam_strain
        # extend mapping table:
        self.mapping_fam["family id"] = strain_fam_index
        self.mapping_fam["original family id"] = strain_fam_labels
        self.mapping_fam["no of strains"] =  np.sum(self.M_fam_strain, axis=1)
        num_members = [x[0].shape[0] for x in self.family_members]
        # see above
        if add_singleton_entries:
            num_members.append(num_of_singletons)
        self.mapping_fam["no of members"] = num_members
        return self.family_members

    def common_strains(self, objects_a, objects_b, filter_no_shared=False):
        """
        Obtain the set of common strains between all pairs from the lists objects_a
        and objects_b.

        The two parameters can be either lists or single instances of the 3 supported
        object types (GCF, Spectrum, MolecularFamily). It's possible to use a single
        object together with a list as well.

        Returns a dict indexed by tuples of (Spectrum/MolecularFamily, GCF), where
        the values are lists of strain indices which appear in both objects, which
        can then be looked up in NPLinker.strains.
        """
        is_list_a = isinstance(objects_a, list)
        is_list_b = isinstance(objects_b, list)

        type_a = type(objects_a[0]) if is_list_a else type(objects_a)
        type_b = type(objects_b[0]) if is_list_b else type(objects_b)

        if type_a == type_b:
            raise Exception('Must supply objects with different types!')
        
        # to keep things slightly simpler, ensure the GCFs are always "b"
        if type_a == GCF:
            type_a, type_b = type_b, type_a
            is_list_a, is_list_b = is_list_b, is_list_a
            objects_a, objects_b = objects_b, objects_a

        if not is_list_a:
            objects_a = [objects_a]
        if not is_list_b:
            objects_b = [objects_b]

        # retrieve object IDs
        ids_b = [gcf.id for gcf in objects_b]
        ids_a = []
        if type_a == Spectrum:
            ids_a = [spec.id for spec in objects_a]
        else:
            # TODO similar treatment to Spectrum/GCF?
            mapping_fam_id = self.mapping_fam["original family id"]
            for fam in objects_a:
                ids_a.append(np.where(mapping_fam_id == int(fam.family_id))[0][0])

        data_a = self.M_spec_strain if type_a == Spectrum else self.M_fam_strain
        data_b = self.M_gcf_strain

        results = {}
        for a, obj_a in enumerate(objects_a):
            for b, obj_b in enumerate(objects_b):
                # just AND both arrays and extract the indices with positive results
                result = np.where(np.logical_and(data_a[ids_a[a]], data_b[ids_b[b]]))[0]
                # if we want to exclude results with no shared strains
                if (filter_no_shared and len(result) > 0) or not filter_no_shared:
                    results[(obj_a, obj_b)] = result

        return results

    def correlation_matrices(self, type='spec-gcf'):
        """
        Collect co-occurrences accros strains:
        IF type='spec-gcf':
            number of co-occurences of spectra and GCFS
            --> Output: M_spec_gcf matrix
        IF type='fam-gcf':
            number of co-occurences of mol.families and GCFS
            --> Output: M_fam_gcf matrix
        """

        # Make selection for scenario spec<->gcf or fam<->gcf
        if type == 'spec-gcf':
            M_type1_strain = self.M_spec_strain
        elif type == 'fam-gcf':
            M_type1_strain = self.M_fam_strain
        elif type == 'spec-bgc' or type == 'fam-bgc':
            raise Exception("Given types are not yet supported... ")
        else:
            raise Exception("Wrong correlation 'type' given. Must be one of 'spec-gcf', 'fam-gcf', ...")

        logger.debug("Calculating correlation matrices of type: {}".format(type))

        # Calculate correlation matrix from co-occurence matrices
        M_type1_gcf, M_type1_notgcf, M_nottype1_gcf, M_nottype1_notgcf = calc_correlation_matrix(M_type1_strain, self.M_gcf_strain)

        # return results:
        if type == 'spec-gcf':
            self.M_spec_gcf = M_type1_gcf
            self.M_spec_notgcf = M_type1_notgcf
            self.M_notspec_gcf = M_nottype1_gcf
            self.M_notspec_notgcf = M_nottype1_notgcf
        elif type == 'fam-gcf':
            self.M_fam_gcf = M_type1_gcf
            self.M_fam_notgcf = M_type1_notgcf
            self.M_notfam_gcf = M_nottype1_gcf
            self.M_notfam_notgcf = M_nottype1_notgcf
        else:
            raise Exception("No correct correlation matrix was created.")

    # class data_links OUTPUT functions
    # TODO add output functions (e.g. to search for mappings of individual specs, gcfs etc.)

class RandomisedDataLinks(DataLinks):
    """Simple subclass of DataLinks that randomly shuffles the strains within
        the strains to spectra matrices (to be used for generating random
        scoring output as point of comparison with real scores) """

    @classmethod
    def from_datalinks(cls, datalinks, find_correlations=True):
        self = cls()
        # check if load_data has been called
        if len(datalinks.M_spec_strain) == 0 or len(datalinks.M_gcf_strain) == 0:
            raise Exception('DataLinks object not initialised (call load_data first)')

        # create copies of the data structures required
        self.M_gcf_strain = datalinks.M_gcf_strain.copy()
        self.M_spec_strain = datalinks.M_spec_strain.copy()
        self.mapping_strain = datalinks.mapping_strain.copy()
        self.mapping_spec = datalinks.mapping_spec.copy()

        # shuffle matrix/matrices
        self._shuffle_cols(self.M_spec_strain)
        # self._shuffle_cols(self.M_gcf_strain)

        # can now run find_correlations with the shuffled matrices
        if find_correlations:
            self.find_correlations()

        return self

    def _shuffle_cols(self, m):
        # shuffle the columns a gcf/spec-strain matrix so that each strain remains 
        # present in the same number of now-random objects.
        # np.random.shuffle only operates on first axis of multi-axis arrays, so
        # take transpose to get columns here instead.
        for col in m.T:
            np.random.shuffle(col)


class LinkLikelihood(object):
    """
    Class to:
    1) create ansd store likelihood matrices (from co-occurences)
    2) select potential calculates for links
    """

    def __init__(self):
        """
        Matrices that store likelihoods of empirically found co-occurence
        Example:
            P_spec_givengcf contains likelihoods P(spec_x|gcf_y),
            which is the probability of finding spec_x given there is a gcf_y
        """

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


    def calculate_likelihoods(self, data_links, type='spec-gcf'):
        """
        Calulate likelihoods from empirically found co-occurences in data
        
        IF type='spec-gcf':
        P(GCF_x | spec_y), P(spec_y | GCF_x),
        P(GCF_x | not spec_y), P(spec_y | not GCF_x)
        IF type='fam-gcf':
        P(GCF_x | fam_y), P(fam_y | GCF_x),
        P(GCF_x | not fam_y), P(fam_y | not GCF_x)
        """

        # Make selection for scenario spec<->gcf or fam<->gcf
        if type == 'spec-gcf':
            M_type1_type2 = data_links.M_spec_gcf
            M_type1_nottype2 = data_links.M_spec_notgcf
            M_nottype1_type2 = data_links.M_notspec_gcf
            M_type1_cond = data_links.M_spec_strain
        elif type == 'fam-gcf':
            M_type1_type2 = data_links.M_fam_gcf
            M_type1_nottype2 = data_links.M_fam_notgcf
            M_nottype1_type2 = data_links.M_notfam_gcf
            M_type1_cond = data_links.M_fam_strain
        elif type == 'spec-bgc' or type == 'fam-bgc':
            raise Exception("Given types are not yet supported... ")
        else:
            raise Exception("Wrong correlation 'type' given. Must be one of 'spec-gcf', 'fam-gcf'...")

        logger.debug("Calculating likelihood matrices of type: {}".format(type))
        # Calculate likelihood matrices using calc_likelihood_matrix()
        P_type2_given_type1, P_type2_not_type1, P_type1_given_type2, \
            P_type1_not_type2 = calc_likelihood_matrix(M_type1_cond,
                                                                              data_links.M_gcf_strain,
                                                                              M_type1_type2,
                                                                              M_type1_nottype2,
                                                                              M_nottype1_type2)
        if type == 'spec-gcf':
            self.P_gcf_given_spec = P_type2_given_type1
            self.P_gcf_not_spec = P_type2_not_type1
            self.P_spec_given_gcf = P_type1_given_type2
            self.P_spec_not_gcf = P_type1_not_type2
        elif type == 'fam-gcf':
            self.P_gcf_given_fam = P_type2_given_type1
            self.P_gcf_not_fam = P_type2_not_type1
            self.P_fam_given_gcf = P_type1_given_type2
            self.P_fam_not_gcf = P_type1_not_type2
        else:
            raise Exception("No correct likelihood matrices were created.")


class LinkFinder(object):
    """
    Class to:
    1) Score potential links based on collected information from:
        DataLinks, LinkLikelihood (and potentially other resources)
    Different scores can be used for this!
    
    2) Rank and output selected candidates
    3) Create output plots and tables
    """

    def __init__(self):
        """
        Create tables of prospective link candidates.
        Separate tables will exist for different linking scenarios, such as
        gcfs <-> spectra OR gcf <-> mol.families
        """
        
        # metcalf scores
        self.metcalf_spec_gcf = []
        self.metcalf_fam_gcf = []

        # likelihood scores
        self.likescores_spec_gcf = []
        self.likescores_fam_gcf = []
        
        # hg scores
        self.hg_spec_gcf = []
        self.hg_fam_gcf = []
        
        # link candidate tables
        self.link_candidates_gcf_spec = []
        self.link_candidates_gcf_fam = []

    def get_scores(self, method, type_):
        if method == 'metcalf':
            if type_ == 'spec-gcf':
                return self.metcalf_spec_gcf
            elif type_ == 'fam-gcf':
                return self.metcalf_fam_gcf
        elif method == 'likescore':
            if type_ == 'spec-gcf':
                return self.likescores_spec_gcf
            elif type_ == 'fam-gcf':
                return self.likescores_fam_gcf
        elif method == 'hg':
            if type_ == 'spec-gcf':
                return self.hg_spec_gcf
            elif type_ == 'fam-gcf':
                return self.hg_fam_gcf

        raise Exception('Unknown method or type (method="{}", type="{}")'.format(method, type_))

    def metcalf_scoring(self, data_links,
                        both=10,
                        type1_not_gcf=-10,
                        gcf_not_type1=0,
                        not_type1_not_gcf=1,
                        type='spec-gcf'):
        """
        Calculate metcalf scores from DataLinks() co-occurence matrices
        """
        
        # TODO expected value scoring stuff
        # Compute the expected values for all possible values of spec and gcf strains
        # we need the total number of strains
        try:
            _,n_strains = self.M_gcf_strain.shape
            expected_metcalf = np.zeros(n_strains+1)
            variance_metcalf = np.zeros(n_strains+1)
            from scipy.stats import hypergeom # maybe needs moving?
            for n in range(n_strains+1):
                for m in range(n_strains+1):
                    max_overlap = min(n,m)
                    min_overlap = max(0,n+m-n_strains) # minimum possible strain overlap
                    expected_value = 0
                    expected_sq = 0
                    for o in range(min_overlap,max_overlap+1):
                        o_prob = hypergeom.pmf(o,n_strains,n,m)
                        # compute metcalf for n strains in type 1 and m in gcf
                        score = o * both
                        score += type1_not_gcf*(n-o)
                        score += gcf_not_type1*(m-o)
                        score += not_type1_not_gcf * (n_strains - (n+m-o))
                        expected_value += o_prob*score
                        expected_sq += o_prob*(score**2)
                    expected_metcalf[n,m] = expected_value
                    variance_metcalf[n,m] = expected_sq - expected_value**2
            # now, we would like an option to take any actual score an subtract the 
            # expected value and then divide by the square root of the variance
            # e.g. if we have a score computed between a type 1 object that has 
            # 3 strains, and a gcf with 6 strains, we should use the expected value 
            # at expected_metcalf[3,6] and sqrt of the variance in the same position
            # 
        except:
            pass

        if type == 'spec-gcf':
            metcalf_scores = np.zeros(data_links.M_spec_gcf.shape)
            metcalf_scores = (data_links.M_spec_gcf * both
                              + data_links.M_spec_notgcf * type1_not_gcf
                              + data_links.M_notspec_gcf * gcf_not_type1
                              + data_links.M_notspec_notgcf * not_type1_not_gcf)
            self.metcalf_spec_gcf = metcalf_scores
            
        elif type == 'fam-gcf':
            metcalf_scores = np.zeros(data_links.M_fam_gcf.shape)
            metcalf_scores = (data_links.M_fam_gcf * both
                              + data_links.M_fam_notgcf * type1_not_gcf
                              + data_links.M_notfam_gcf * gcf_not_type1
                              + data_links.M_notfam_notgcf * not_type1_not_gcf)
            
            self.metcalf_fam_gcf = metcalf_scores
        return metcalf_scores

    def hg_scoring(self, data_links, type='spec-gcf'):
        """
        Calculate metcalf scores from DataLinks() co-occurence matrices
        """

        # NOTE:can't use the correlation matrices directly for this scoring method because
        # it seems to require more inclusive counts of the strains in each object. 
        # Instead of "number of strains only in GCF", it requires "number of strains in the
        # GCF PLUS the number shared between the GCF and the other object". 
        # e.g. if a spectrum has 3 strains, a GCF has 1 strain and there is 1 shared strain, 
        # M_spec_gcf will correctly contain "1", but M_type1_notgcf will contain "2" instead
        # of "3", because the spectrum only has 2 distinct strains vs the GCF.
        # To fix this the M_spec_gcf/M_fam_gcf matrix can just be added onto the others to give
        # the correct totals.


        if type == 'spec-gcf':
            num_strains = np.ones(data_links.M_spec_gcf.shape) * data_links.M_gcf_strain.shape[1]
            overlap_counts = data_links.M_spec_gcf
            gcf_counts = overlap_counts + data_links.M_notspec_gcf
            spec_counts = overlap_counts + data_links.M_spec_notgcf
            hg_scores = hypergeom.sf(overlap_counts, num_strains, gcf_counts, spec_counts, loc=1)
            self.hg_spec_gcf = hg_scores
        elif type == 'fam-gcf':
            num_strains = np.ones(data_links.M_fam_gcf.shape) * data_links.M_gcf_strain.shape[1]
            overlap_counts = data_links.M_fam_gcf
            gcf_counts = overlap_counts + data_links.M_notfam_gcf
            fam_counts = overlap_counts + data_links.M_fam_notgcf
            hg_scores = hypergeom.sf(overlap_counts, num_strains, gcf_counts, fam_counts, loc=1)
            self.hg_fam_gcf = hg_scores

        return hg_scores

    def likelihood_scoring(self, data_links, likelihoods,
                        alpha_weighing=0.5,
                        type='spec-gcf'):
        """
        Calculate likelihood scores from DataLinks() co-occurence matrices.
        
        Idea:
            Score reflect the directionality BGC-->compound-->spectrum, which
            suggests that the most relevant likelihoods are:
            P(gcf|type1) - If type1 is the result of only one particular gene cluster,
                            this value should be high (close or equal to 1)
            P(type1|not gcf) - Following the same logic, this value should be very
                                small or 0 ("no gene cluster, no compound")
                                
        Score:
            Score = P(gcf|type1) * (1 - P(type1|not gcf) * weighing function
            
            weighing function is here a function of the number of strains they co-occur.
            weighing function = (1 - exp(-alpha_weighing * num_of_co-occurrences)
        """
        
        if type == 'spec-gcf':
            likelihood_scores = np.zeros(data_links.M_spec_gcf.shape)
            likelihood_scores = (likelihoods.P_gcf_given_spec
                                 * (1 - likelihoods.P_spec_not_gcf)
                                 * (1 - np.exp(-alpha_weighing * data_links.M_spec_gcf))
                                 )
            
            self.likescores_spec_gcf = likelihood_scores
            
        elif type == 'fam-gcf':
            likelihood_scores = np.zeros(data_links.M_fam_gcf.shape)
            likelihood_scores = (likelihoods.P_gcf_given_fam
                                 * (1 - likelihoods.P_fam_not_gcf)
                                 * (1 - np.exp(-alpha_weighing * data_links.M_fam_gcf))
                                 )
            
            self.likescores_fam_gcf = likelihood_scores
        return likelihood_scores


    def select_link_candidates(self, data_links, likelihoods,
                               P_cutoff=0.8,
                               main_score='likescore',
                               score_cutoff=0,
                               type='fam-gcf'):
        """
        Look for potential best candidate for links between
        IF type='spec-gcf': GCFs and spectra
        IF type='fam-gcf': GCFs and mol.families
        
        Parameters
        ----------
        data_links: DataLinks() object
            containing co-occurence matrices
        likelihood: LinkLikelihood() object
            containing co-occurence likelihoods
        P_cutoff: float
            Thresholds to conly consider candidates for which:
            P_gcf_given_type1 >= P_cutoff
        main_score: str
            Which main score to use ('metcalf', 'likescore')
        score_cutoff:
            Thresholds to conly consider candidates for which:
            score >= score_cutoff
        """

        # Select scenario: spec<->gcf or fam<->gcf
        if type == 'spec-gcf':
            P_gcf_given_type1 = likelihoods.P_gcf_given_spec
            P_gcf_not_type1 = likelihoods.P_gcf_not_spec
            P_type1_given_gcf = likelihoods.P_spec_given_gcf
            P_type1_not_gcf = likelihoods.P_spec_not_gcf
            M_type1_gcf = data_links.M_spec_gcf
            metcalf_scores = self.metcalf_spec_gcf
            likescores = self.likescores_spec_gcf
            index_names = ["spectrum_id", "GCF id", "P(gcf|spec)", "P(spec|gcf)",
                             "P(gcf|not spec)", "P(spec|not gcf)",
                             "co-occur in # strains",
                             "metcalf score", "likelihood score",
                             "HG prob", "link prob", "link prob specific"]

        elif type == 'fam-gcf':
            P_gcf_given_type1 = likelihoods.P_gcf_given_fam
            P_gcf_not_type1 = likelihoods.P_gcf_not_fam
            P_type1_given_gcf = likelihoods.P_fam_given_gcf
            P_type1_not_gcf = likelihoods.P_fam_not_gcf
            M_type1_gcf = data_links.M_fam_gcf
            metcalf_scores = self.metcalf_fam_gcf
            likescores = self.likescores_fam_gcf
            index_names = ["family_id", "GCF id", "P(gcf|fam)", "P(fam|gcf)",
                             "P(gcf|not fam)", "P(fam|not gcf)",
                             "co-occur in # strains",
                             "metcalf score", "likelihood score",
                             "HG prob", "link prob", "link prob specific"]

        elif type == 'spec-bgc' or type == 'fam-bgc':
            raise Exception("Given types are not yet supported... ")
        else:
            raise Exception("Wrong correlation 'type' given. Must be one of 'spec-gcf', 'fam-gcf'...")

        dim1, dim2 = P_gcf_given_type1.shape

        # PRE-SELECTION: 
        # Select candidates with P_gcf_given_spec >= P_cutoff AND score >= score_cutoff        
        if main_score == 'likescore':
            candidate_ids = np.where((P_gcf_given_type1[:, :] >= P_cutoff) & (likescores >= score_cutoff))
        elif main_score == 'metcalf':
            candidate_ids = np.where((P_gcf_given_type1[:, :] >= P_cutoff) & (metcalf_scores >= score_cutoff))
        else:
            raise Exception('Wrong scoring type given. Must be one of: {}'.format(SCORING_METHODS))
            
        logger.debug(candidate_ids[0].shape[0], " candidates selected with ",
              index_names[2], " >= ", P_cutoff, " and a link score >= ", score_cutoff, ".")
        
        link_candidates = np.zeros((12, candidate_ids[0].shape[0]))
        link_candidates[0, :] = candidate_ids[0]  # spectrum/fam number
        link_candidates[1, :] = candidate_ids[1]  # gcf id
        link_candidates[2, :] = P_gcf_given_type1[candidate_ids]
        link_candidates[3, :] = P_type1_given_gcf[candidate_ids]
        link_candidates[4, :] = P_gcf_not_type1[candidate_ids]
        link_candidates[5, :] = P_type1_not_gcf[candidate_ids]
        link_candidates[6, :] = M_type1_gcf[candidate_ids]
        link_candidates[7, :] = metcalf_scores[candidate_ids]
        link_candidates[8, :] = likescores[candidate_ids]
        
        # Calculate probability to find similar link by chance
        Nx_list = data_links.mapping_gcf["no of strains"]
        if type == 'spec-gcf':
            Ny_list = data_links.mapping_spec["no of strains"]
        elif type == 'fam-gcf':
            Ny_list = data_links.mapping_fam["no of strains"]
        
        # Calculate probabilities of finding a spectrum in a certain strain
        P_str = np.array(data_links.mapping_strain["no of spectra"])
        P_str = P_str/np.sum(P_str)
            
        num_strains = data_links.M_gcf_strain.shape[1]
        
        # Calculate the hypergeometric probability (as before)
        for i in range(link_candidates.shape[1]):
             link_candidates[9, i] = pair_prob_hg(link_candidates[6,i],
                                                            num_strains,
                                                            Nx_list[link_candidates[1,i]],
                                                            Ny_list[int(link_candidates[0,i])])

        # Calculate the GCF specific probability
        for i in range(link_candidates.shape[1]):
            id_spec = int(link_candidates[0,i])
            id_gcf = int(link_candidates[1,i])
            
            # find set of strains which contain GCF with id link_candidates[1,i] 
            XG = np.where(data_links.M_gcf_strain[id_gcf , :] == 1)[0]
                                                           
            link_candidates[10,i] = pair_prob_approx(P_str, XG,
                                                            int(Ny_list[id_spec]),
                                                            int(link_candidates[6,i]))
            # Calculate the link specific probability
            # Find strains where GCF and spectra/family co-occur
            if type == 'spec-gcf':
                XGS = np.where((data_links.M_gcf_strain[id_gcf, :] == 1) & (data_links.M_spec_strain[id_spec, :] == 1))[0]
            elif type == 'fam-gcf':
                XGS = np.where((data_links.M_gcf_strain[id_gcf, :] == 1) & (data_links.M_fam_strain[id_spec, :] == 1))[0]
            link_candidates[11,i] = link_prob(P_str, XGS,
                                                           int(Nx_list[id_gcf]),
                                                           int(Ny_list[id_spec]), num_strains)
            
        # Transform into pandas Dataframe (to avoid index confusions):
        link_candidates_pd = pd.DataFrame(link_candidates.transpose(1,0), columns = index_names)
        
        # add other potentially relevant knowdledge
        # If this will grow to more collected information -> create separate function/class
        bgc_class = []
        for i in link_candidates_pd["GCF id"].astype(int):
            bgc_class.append(data_links.mapping_gcf["bgc class"][i])
        link_candidates_pd["BGC class"] = bgc_class

        # Change some columns to int
        link_candidates_pd.iloc[:,[0,1,6,7]] = link_candidates_pd.iloc[:,[0,1,6,7]].astype(int)
        
        # return results
        if type == 'spec-gcf':
            self.link_candidates_gcf_spec = link_candidates_pd
        elif type == 'fam-gcf':
            self.link_candidates_gcf_fam = link_candidates_pd
        else:
            raise Exception("No candidate selection was created.")
        return link_candidates_pd
    
    
    def get_links(self, data_links, input_object,
                  main_score='likescore',
                  score_cutoff=0.5):
        """
        Output likely links for 'input_object'
        
        Parameters
        ----------
        input_objects: object()
            Object or list of objects of either class: spectra, families, or GCFs
        main_score: str
            Which main score to use ('metcalf', 'likescore')
        score_cutoff:
            Thresholds to conly consider candidates for which:
            score >= score_cutoff
        """

        if main_score not in SCORING_METHODS:
            raise Exception('Wrong scoring type given. Must be one of: {}'.format(SCORING_METHODS))

        # Check if input is list:
        if isinstance(input_object, list):
            query_size = len(input_object)
        else:
            input_object = [input_object]
            query_size = 1

        # Check type of input_object:
        # If GCF:
        if isinstance(input_object[0], GCF):
            input_type = "gcf"
            link_levels = [0, 1]
            
            # Get necessary ids
            input_ids = np.zeros(query_size)
            for i, gcf in enumerate(input_object):
                input_ids[i] = gcf.id
            input_ids =  input_ids.astype(int)
            
            if main_score == 'likescore':
                likescores = [self.likescores_spec_gcf[:, input_ids],
                          self.likescores_fam_gcf[:, input_ids]]
            elif main_score == 'metcalf':
                metcalf_scores = [self.metcalf_spec_gcf[:, input_ids],
                              self.metcalf_fam_gcf[:, input_ids]]
            elif main_score == 'hg':
                hg_scores = [self.hg_spec_gcf[:, input_ids],
                         self.hg_fam_gcf[:, input_ids]]

        # If Spectrum:
        elif isinstance(input_object[0], Spectrum):
            
            # Get necessary ids
            input_ids = np.zeros(query_size)
            mapping_spec_id = data_links.mapping_spec["original spec-id"]
            for i, spectrum in enumerate(input_object):
                input_ids[i] = np.where(mapping_spec_id == int(spectrum.id))[0]
            input_ids =  input_ids.astype(int)

            input_type = "spec"
            link_levels = [0]
            if main_score == 'likescore':
                likescores = [self.likescores_spec_gcf[input_ids, :],
                              []]
            elif main_score == 'metcalf':
                metcalf_scores = [self.metcalf_spec_gcf[input_ids, :],
                                  []]
            elif main_score == 'hg':
                hg_scores = [self.hg_spec_gcf[input_ids, :],
                             []]
        # If MolecularFamily:                     
        elif isinstance(input_object[0], MolecularFamily):
            
            # Get necessary ids
            # TODO: include Singletons, maybe optinal
            input_ids = np.zeros(query_size)
            mapping_fam_id = data_links.mapping_fam["original family id"]
            for i, family in enumerate(input_object):
                input_ids[i] = np.where(mapping_fam_id == int(family.family_id))[0]
            input_ids =  input_ids.astype(int)
            
            input_type = "fam"
            link_levels = [1]
            if main_score == 'likescore':
                likescores = [[],
                              self.likescores_fam_gcf[input_ids, :]]
            elif main_score == 'metcalf':
                metcalf_scores = [[],
                                  self.metcalf_fam_gcf[input_ids, :]]
            elif main_score == 'hg':
                hg_scores = [[],
                             self.hg_fam_gcf[input_ids, :]]
        else:
            raise Exception("Input_object must be Spectrum, MolecularFamily, or GCF object (single or list).")

        links = []
        for linklevel in link_levels:
            if score_cutoff is not None:
                if main_score == 'likescore':
                    candidate_ids = np.where(likescores[linklevel] >= score_cutoff)
                elif main_score == 'metcalf':
                    candidate_ids = np.where(metcalf_scores[linklevel] >= score_cutoff)
                elif main_score == 'hg':
                    candidate_ids = np.where(hg_scores[linklevel] >= score_cutoff)
                else:
                    # should never happen
                    raise Exception('Unknown scoring type! "{}"'.format(main_score))
            else:
                # TODO is this best way to get same output as above code?
                # to keep the remainder of the method identical in the case of no cutoff
                # being supplied, while still returning all the candidate links, I'm
                # currently abusing np.where like this
                candidate_ids = np.where(metcalf_scores[linklevel] != np.nan)

            link_candidates = np.zeros((3, candidate_ids[0].shape[0]))
            
            if query_size == 1:
                link_candidates[0, :] = input_ids  # input id
                if input_type == 'gcf':
                    link_candidates[1, :] = candidate_ids[0].astype(int)  # output id
                else:
                    link_candidates[1, :] = candidate_ids[1].astype(int)  # output id
            else:
                if input_type == 'gcf':
                    link_candidates[0, :] = input_ids[candidate_ids[1]]  # input id
                    link_candidates[1, :] = candidate_ids[0].astype(int)  # output id
                else:
                    link_candidates[0, :] = input_ids[candidate_ids[0]]  # input id
                    link_candidates[1, :] = candidate_ids[1].astype(int)  # output id
            
            if main_score == 'likescore':
                link_candidates[2, :] = likescores[linklevel][candidate_ids]
            elif main_score == 'metcalf':
                link_candidates[2, :] = metcalf_scores[linklevel][candidate_ids]
            elif main_score == 'hg':
                link_candidates[2, :] = hg_scores[linklevel][candidate_ids]

            links.append(link_candidates)

        return links




    def create_cytoscape_files(self, data_links,
                               network_filename,
                               link_type='fam-gcf',
                               score_type='metcalf'):
        """
        Create network file for import into Cytoscape.
        Network file will be generated using networkx.
        The type of network created here is a bipartite network.
            mass spec side --> bipartite = 0
            gene cluster side --> bipartite = 1
        Output format is a graphml file.
        
        Parameters
        ----------
        data_links: DataLinks() object
            containing co-occurence matrices
        likelihood: LinkLikelihood() object
            containing co-occurence likelihoods
        network_filename: str
            Filename to save generated model as graphml file.
        """
        
        import networkx as nx
        NPlinker_net = nx.Graph()
        
        if link_type == 'fam-gcf':
            link_candidates = self.link_candidates_gcf_fam
            type1str = 'family_id'
        elif link_type == 'spec-gcf':
            link_candidates = self.link_candidates_gcf_spec
            type1str = 'spectrum_id'
        else:
            raise Exception("Wrong link-type given.")
            
        # Select score type
        if score_type == 'metcalf':
            scorestr = 'metcalf score'
        elif score_type == 'likescore':
            scorestr = 'likelihood score'
        else:
            raise Exception("Wrong score_type given. Must be one of: 'metcalf', 'likescore' .")
        
        # Add nodes (all partners from link_candidate table):
        # mass spec side --> bipartite = 0
        type1_names = []
        for type1 in link_candidates[type1str].astype(int):
            type1_names.append(type1str+str(type1))
            
        type1_names_unique = list(set(type1_names))
        NPlinker_net.add_nodes_from(type1_names_unique, bipartite=0)
        
        # gene cluster side --> bipartite = 1
        type2_names = []
        for type2 in link_candidates['GCF id']:
            type2_names.append("GCF_"+str(type2))
     
        type2_names_unique = list(set(type2_names))
        NPlinker_net.add_nodes_from(type2_names_unique, bipartite=1)
        
        # Add edges:
        for i in range(0, link_candidates.shape[0]):
            NPlinker_net.add_edge(type1_names[i],
                                  type2_names[i],
                                  weight=float(link_candidates[scorestr][i]))
            
        # Add edges between molecular family members
        if link_type == 'spec-gcf':
            type1_unique = (np.unique(link_candidates[type1str])).astype(int)
            map_spec_fam = data_links.mapping_spec["fam-id"]
            for type1 in type1_unique:
                if map_spec_fam[type1] > 0:  # if no singleton
                    members = data_links.family_members[int(map_spec_fam[type1])][0]
                    
                    # select other family members if among link candidates
                    members_present = [x for x in list(members) if x in list(type1_unique)]
                    for member in members_present:
                        NPlinker_net.add_edge(type1str+str(type1), type1str+str(member))
                
        # export graph for drawing (e.g. using Cytoscape)
        nx.write_graphml(NPlinker_net, network_filename)


    def plot_candidates(self, P_cutoff=0.8, 
                        score_type='likescore', 
                        score_cutoff=0, 
                        type='fam-gcf'):
        """
        Plot best rated correlations between gcfs and spectra/families
        plot in form of seaborn clustermap
        """
        
        # Select score type
        if score_type == 'metcalf':
            scorestr = 'metcalf score'
        elif score_type == 'likescore':
            scorestr = 'likelihood score'
        else:
            raise Exception("Wrong score_type given. Must be one of: 'metcalf', 'likescore' .")
        
        if type == 'spec-gcf':
            link_candidates = self.link_candidates_gcf_spec
            selected_ids = np.where((link_candidates["P(gcf|spec)"] > P_cutoff) & 
                        (link_candidates[scorestr] > score_cutoff))[0]
        elif type == 'fam-gcf':   
            link_candidates = self.link_candidates_gcf_fam
            selected_ids = np.where((link_candidates["P(gcf|fam)"] > P_cutoff) & 
                        (link_candidates[scorestr] > score_cutoff))[0]
        else:
            raise Exception("Wrong correlation 'type' given.")   
    
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
            M_links[x, y] = link_candidates[scorestr][i] 
    
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
                               vmin=0, vmax=np.max(link_candidates[scorestr]), 
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

