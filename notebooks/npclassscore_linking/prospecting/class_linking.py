#!/usr/bin/env python3
'''
Initial code for NPClassScore
'''
import sys
import os
import glob
import pandas as pd
from collections import defaultdict, Counter

sys.path.append('../../src')
from nplinker.nplinker import BGC, GCF, Spectrum
from nplinker.nplinker import NPLinker


class Class_links:
    '''Holds all info concerning class links (based on known bgc-compound links in MIBiG)
    '''
    def __init__(self, mibig_classes_file):
        self._mibig_classes_file = mibig_classes_file
        self._read_mibig_classes()
        self._get_class_counts()
        self._get_scoring_tables()
        pd.options.display.float_format = "{:,.3f}".format  # adjust pd formatting

        self._bigscape_mibig_conversion = {
            'PKSI': 'Polyketide', 'PKSother': 'Polyketide',
            'NRPS': 'NRP', 'RiPPs': 'RiPP', 'Saccharides': 'Saccharide',
            'Others': 'Other', 'Terpene': 'Terpene', 'PKS-NRP_Hybrids': 'PKS-NRP_Hybrids'}

        self._as_conversion = {
            'NAGGN': 'other', 'NAPAA': 'other', 'RRE-containing': 'bacteriocin',
            'RiPP-like': 'bacteriocin', 'cf_fatty_acid': "fatty_acid", 'cf_putative': 'other',
            'cf_saccharide': 'saccharide', 'guanidinotides': 'fused', 'lanthipeptide-class-i': 'lanthipeptide',
            'lanthipeptide-class-ii': 'lanthipeptide', 'lanthipeptide-class-iii': 'lanthipeptide',
            'lanthipeptide-class-iv': 'lanthipeptide', 'lanthipeptide-class-v': 'lanthipeptide',
            'lantipeptide': 'lanthipeptide', 'linaridin': 'lanthipeptide', 'lipolanthine': 'lanthipeptide',
            'nrps': 'NRPS', 'otherks': 'hglE-KS', 'prodigiosin': 'other', 'pyrrolidine': 'other',
            'ranthipeptide': 'bacteriocin', 'redox-cofactor': 'other', 't1pks': 'T1PKS', 't2pks': 'T2PKS',
            't3pks': 'T3PKS', 'thioamide-NRP': 'other', 'thioamitides': 'bacteriocin',
            'transatpks': 'transAT-PKS'}

    def get_gcf_as_classes(self, gcf, cutoff = 0.5):
        """Get antismash classes for a gcf if antismash class occurs in more than <cutoff> of gcf

        Args:
            - gcf: GCF NPLinker object
            - cutoff: float - fraction of the GCF that needs to contain the class for it to be counted
        Returns:
            List of str, the antismash classes present in the gcf
        """
        # todo: move to GCF object?
        gcf_size = len(gcf.bgcs)
        unlist_all_products = [product for bgc in gcf.bgcs for product in bgc.product_prediction.split('.')]
        sorted_as_classes = Counter(unlist_all_products).most_common()
        # keep if in more than half of bgcs?
        cutoff = 0.5
        size_cutoff = gcf_size * cutoff
        filtered_as_classes = []
        for product in sorted_as_classes:
            if product[1] >= size_cutoff:
                filtered_as_classes.append(product[0])
        return filtered_as_classes

    def convert_as_classes(self, init_as_classes: list):
        """Convert AS classes to class names that are in scoring table

        Args:
            - init_as_classes: list of str, the initial antismash class names
        Returns:
            List of str: converted antismash classes with _as_conversion_table
        """
        as_classes = []
        for as_class in init_as_classes:
            as_conversion = self.as_conversion.get(as_class)
            if as_conversion:
                as_classes.append(as_conversion)
            else:
                as_classes.append(as_class)
        return as_classes

    def _read_mibig_classes(self):
        # read mibig file to dict of list {chem_id: [bgc_classes, chem_classes]}
        classes_dict = {}
        with open(self._mibig_classes_file) as inf:
            header = inf.readline()
            for line in inf:
                elems = line.strip().split("\t")
                chem_id = elems.pop(0)
                class_base = elems.pop(0).split(',')
                classes = [cls.partition(':')[0] for cls in class_base]
                sub_classes = [cls for cls in class_base if cls.split(":")[1]]
                as_classes = elems.pop(0).split(',')

                bgc_classes = [classes, sub_classes, as_classes]
                chem_classes = [chem_cls.split('; ') for chem_cls in elems[2:]]
                classes_dict[chem_id] = [bgc_classes, chem_classes]
        self._mibig_classes = classes_dict
        # add header info
        s_h = header.strip().split('\t')

        self._bgc_class_names = ['mibig_classes']+s_h[1:3]
        self._chem_class_names = s_h[5:]

        return self._mibig_classes

    def _get_class_counts(self):
        # aggregate pairwise class matrices for all compounds

        def _rec_dd():
            """Initialises a recurring defaultdict"""
            return defaultdict(_rec_dd)

        result = _rec_dd()
        # loop through each mibig compound
        for mibig_chem_id, (bgc_classes, chem_classes) in self._mibig_classes.items():
        # get all combinations of classes for this compound
            for i, bgc_cat in enumerate(self.bgc_class_names):
                init_bgc_class = bgc_classes[i]
                if not init_bgc_class or init_bgc_class == ['']:
                    continue

                bgc_class = init_bgc_class[:]  # if no exceptions, just assign classes

                # do some cleanup for mibig classes
                if bgc_cat == "mibig_classes":
                    # group pks-nrp hybrids for MIBiG classes
                    hyb_count = len([1 for init_bgc_c in init_bgc_class \
                                     if any([test in init_bgc_c.lower() for test in ['nrp', 'pks', 'polyketide']])])
                    if hyb_count >= 2:
                        # if hybrid, reconstruct the bgc_class
                        bgc_class = []
                        bgc_class.append("PKS-NRP_Hybrids")
                        # append other classes if there are more
                        for init_bgc_c in init_bgc_class:
                            if not any([test in init_bgc_c.lower() for test in ['nrp', 'pks', 'polyketide']]):
                                bgc_class.append(init_bgc_c)

                    # replace Alkaloid with Other in bgc_class
                    bgc_class = ["Other" if bgc_c == "Alkaloid" else bgc_c for bgc_c in bgc_class]

                for j, chem_cat in enumerate(self.chem_class_names):
                    chem_class = chem_classes[j]
                    if not chem_class or chem_class == ['']:
                        continue

                    for bgc_c in bgc_class:
                        for chem_c in chem_class:
                            try:
                                result[bgc_cat][chem_cat][bgc_c][chem_c] += 1
                            except TypeError:
                                result[bgc_cat][chem_cat][bgc_c][chem_c] = 1
        self._class_count_dict = result
        return result

    def _get_scoring_tables(self):
        # makes dict of pd.DataFrames
        # read resulting tables column to row
        # bgc -> chem: d[bgc_key][chem_key][bgc_class][chem_class]
        # and vice versa for chem -> bgc
        class_linking_tables = {}
        class_linking_counts = {}  # store the counts in df/get rid of defaultdicts
        for bgc_key, bgc_chem_counts in self._class_count_dict.items():
            for chem_key, counts in bgc_chem_counts.items():
                # init entries in dict
                if not bgc_key in class_linking_tables:
                    class_linking_tables[bgc_key] = {}
                    class_linking_counts[bgc_key] = {}
                if not chem_key in class_linking_tables:
                    class_linking_tables[chem_key] = {}
                    class_linking_counts[chem_key] = {}
                # add linking tables as DataFrames
                counts_df = pd.DataFrame.from_dict(counts, dtype=int)
                class_linking_tables[bgc_key][chem_key] = (counts_df/counts_df.sum(axis=0)).fillna(0)
                class_linking_counts[bgc_key][chem_key] = counts_df.fillna(0)
                class_linking_tables[chem_key][bgc_key] = (counts_df.T/counts_df.sum(axis=1)).fillna(0)
                class_linking_counts[chem_key][bgc_key] = counts_df.T.fillna(0)
        self._class_links = class_linking_tables
        self._class_links_counts = class_linking_counts
        return class_linking_tables

    @property
    def class_links(self):
        return self._class_links

    @property
    def class_links_counts(self):
        return self._class_links_counts

    @property
    def bgc_class_names(self):
        return self._bgc_class_names

    @property
    def chem_class_names(self):
        return self._chem_class_names

    @property
    def bigscape_mibig_conversion(self):
        return self._bigscape_mibig_conversion

    @property
    def as_conversion(self):
        return self._as_conversion


class Canopus_results:
    """Class for storing canopus results

    The two input files from input_dir are read for the spectra and molfams, respectively:
        -cluster_index_classifications.txt
        -component_index_classifications.txt
    """
    def __init__(self, root_dir):
        """Read the class info from root_dir/canopus

        Args:
            root_dir: str, root_dir of nplinker project
        """
        spectra_classes_names, spectra_classes = self._read_spectra_classes(root_dir)
        self._spectra_classes = spectra_classes
        self._spectra_classes_names = spectra_classes_names
        self._spectra_classes_names_inds = {elem:i for i,elem in enumerate(spectra_classes_names)}

        molfam_classes_names, molfam_classes = self._read_molfam_classes(root_dir)
        self._molfam_classes = molfam_classes
        self._molfam_classes_names = molfam_classes_names
        self._molfam_classes_names_inds = {elem:i for i,elem in enumerate(molfam_classes_names)}

    def _read_spectra_classes(self, root_dir):
        """Read canopus classes for spectra, return classes_names, classes

        Args:
            root_dir: str, root_dir of nplinker project
        Returns:
            Tuple of:
            - ci_classes_names: list of str - the names of each different level
            - ci_classes: dict of {str, lists of tuple(str, float)} - per spectrum (key) the classes for each level
                where each level is a list of (class_name, score) sorted on best choice so index 0 is the best
                class prediction for a level
        """
        input_file = os.path.join(
            root_dir, 'canopus', 'cluster_index_classifications.txt')

        ci_classes = {}  # make a dict {ci: [[(class,score)]]}
        ci_classes_header = None
        ci_classes_names = []

        if os.path.isfile(input_file):
            with open(input_file) as inf:
                ci_classes_header = inf.readline().strip().split("\t")
                for line in inf:
                    line = line.strip('\n').split("\t")
                    classes_list = []
                    for lvl in line[3:]:
                        lvl_list = []
                        for l_class in lvl.split("; "):
                            if l_class:
                                l_class = l_class.split(":")
                                c_tup = tuple([l_class[0], float(l_class[1])])
                            else:
                                c_tup = None  # default value for class value
                            lvl_list.append(c_tup)
                        classes_list.append(lvl_list)
                    ci_classes[line[1]] = classes_list
        if ci_classes_header:
            #  todo: rename the output from the canopus script directly
            ci_classes_names = [f"cf_{elem}" for elem in ci_classes_header[3:-3]] +\
                   [f"npc_{elem}" for elem in ci_classes_header[-3:]]
        return ci_classes_names, ci_classes

    def _read_molfam_classes(self, root_dir):
        """Read canopus classes for molfams, return classes_names, classes

        Args:
            root_dir: str, root_dir of nplinker project
        Returns:
            Tuple of:
            - compi_classes_names: list of str - the names of each different level
            - compi_classes: dict of {str: lists of tuple(str, float)} - per molfam (key) the classes for each level
                where each level is a list of (class_name, fraction) sorted on best choice so index 0 is the best
                class prediction for a level
        """
        input_file = os.path.join(
            root_dir, 'canopus', 'component_index_classifications.txt')

        compi_classes = {}  # make a dict {compi: [[(class,score)]]}
        compi_classes_header = None
        compi_classes_names = []

        if os.path.isfile(input_file):
            with open(input_file) as inf:
                compi_classes_header = inf.readline().strip().split("\t")
                for line in inf:
                    line = line.strip('\n').split("\t")
                    classes_list = []
                    for lvl in line[2:]:
                        lvl_list = []
                        for l_class in lvl.split("; "):
                            if l_class:
                                l_class = l_class.split(":")
                                c_tup = tuple([l_class[0], float(l_class[1])])
                            else:
                                c_tup = None  # default value for class value
                            lvl_list.append(c_tup)
                        classes_list.append(lvl_list)
                    compi_classes[line[0]] = classes_list
        if compi_classes_header:
            #  todo: rename the output from the canopus script directly
            compi_classes_names = [f"cf_{elem}" for elem in compi_classes_header[2:-3]] +\
                      [f"npc_{elem}" for elem in compi_classes_header[-3:]]
        return compi_classes_names, compi_classes

    @property
    def spectra_classes(self):
        return self._spectra_classes

    @property
    def spectra_classes_names(self):
        return self._spectra_classes_names

    @property
    def spectra_classes_names_inds(self):
        return self._spectra_classes_names_inds

    @property
    def molfam_classes(self):
        return self._molfam_classes

    @property
    def molfam_classes_names(self):
        return self._molfam_classes_names

    @property
    def molfam_classes_names_inds(self):
        return self._molfam_classes_names_inds


class MolNetEnhancer_results:
    """Class for storing MolNetEnhancer results

    The input file for ClassyFire results is read from the molnetenhancer directory:
        - ClassyFireResults_Network.txt
    """
    def __init__(self, root_dir):
        """Read the class info from file in root_dir/molnetenhancer/

        Args:
            root_dir: str, root_dir of nplinker project
        """
        cf_classes_names, molfam_classes, spectra2molfam = self._read_cf_classes(root_dir)
        self._spectra2molfam = spectra2molfam
        self._molfam_classes = molfam_classes
        self._spectra_classes_names = cf_classes_names  # if NPC gets implemented, add here
        self._spectra_classes_names_inds = {elem:i for i,elem in enumerate(cf_classes_names)}

    def _read_cf_classes(self, root_dir):
        r"""Read ClassyFireResults_Network.txt in molnetenhancer dir

        Args:
            root_dir: str, root_dir of nplinker project
        Returns:
            tuple of:
            -list of str - names of the classes in order
            -dict of {str: [(str, float)]} - linking molfams to (classes, scores) in order of names,
                singleton families are denoted with S[\d]+
            -dict of {str:str} - linking spectra to molfams
        """
        columns = []
        mne_component_dict = {}
        mne_cluster2component = {}
        input_file = os.path.join(root_dir, "molnetenhancer", "ClassyFireResults_Network.txt")
        if not os.path.isfile(input_file):
            print("no molnetenhancer input")  # todo: add to logger
            return columns, mne_component_dict, mne_cluster2component
        with open(input_file) as inf:
            header = inf.readline().strip().split("\t")
            columns = ['cf_direct_parent' if col == 'CF_Dparent' else col.lower()
                       for i, col in enumerate(header[3:]) if i%2 == 0]
            for line in inf:
                line = line.strip('\n').split("\t")
                cluster = line.pop(0)
                component = line.pop(0)
                nr_nodes = line.pop(0)
                class_info = []
                for i in range(0, len(line), 2):
                    class_tup = (line[i], float(line[i+1]))
                    class_info.append(class_tup)
                if not component in mne_component_dict:
                    mne_component_dict[component] = class_info
                mne_cluster2component[cluster] = component

        return columns, mne_component_dict, mne_cluster2component

    def spectra_classes(self, spectrum_id):
        """Return classes by relating spectrum_id in the molfam_classes

        Args:
            spectrum_id: int/str, spectrum_id - ints will be converted to str
        """
        classes = []
        if isinstance(spectrum_id, int):
            spectrum_id = str(spectrum_id)
        molfam_id = self.spectra2molfam.get(spectrum_id)
        if molfam_id:
            classes = self.molfam_classes.get(molfam_id)
        return classes

    @property
    def spectra2molfam(self):
        return self._spectra2molfam

    @property
    def spectra_classes_names(self):
        return self._spectra_classes_names

    @property
    def spectra_classes_names_inds(self):
        return self._spectra_classes_names_inds

    @property
    def molfam_classes(self):
        return self._molfam_classes


# dummy class before integrating it into the NPLinker class
class NPLinker_classes(NPLinker):
    def __init__(self, userconfig=None):
        super().__init__(userconfig)

    def read_class_info(self):
        """Should include this in load_data/loader()

        returns a Canopus and MNE object with info about classes, and the Class links object
        """
        mibig_classes_file = glob.glob(os.path.join(self.root_dir, "MIBiG*_compounds_with_AS_BGC_CF_NPC_classes.txt"))[0]
        self._class_links = Class_links(mibig_classes_file)

        class_predict_options = []
        self._canopus = Canopus_results(self.root_dir)
        if self._canopus.spectra_classes:
            print("CANOPUS classes loaded succesfully")
            class_predict_options.append('canopus')
        self._molnetenhancer = MolNetEnhancer_results(self.root_dir)
        if self._molnetenhancer.spectra2molfam:
            print("MolNetEnhancer classes loaded succesfully")
            class_predict_options.append('molnetenhancer')
        if class_predict_options:
            class_predict_options = ['mix', 'main'] + class_predict_options
        self._class_predict_options = class_predict_options

    def class_linking_score(self, obj, target):
        """Return sorted class link scores for scoring obj and target

        The input objects can be any mix of the following NPLinker types:
            - BGC
            - GCF
            - Spectrum
            - MolecularFamily

        Args:
            - obj: one of the possible input objects
            - target: one of the possible input objects
        Returns:
            List of tuple of
            (score, obj_class_lvl, target_class_lvl, obj_class, target_class)
            (float, str, str, str, str)
        """
        # assess what is obj and target
        spec_like = obj
        bgc_like = target
        if isinstance(obj, BGC) or isinstance(obj, GCF):
            bgc_like = obj
            spec_like = target

        # assess if bgc or gcf and if spectrum or molfam
        is_spectrum = False
        is_bgc = False
        if isinstance(spec_like, Spectrum):
            is_spectrum = True
        if isinstance(bgc_like, BGC):
            is_bgc = True

        # gather correct classes based on input, dict for bgcs and list for spec
        if is_bgc:
            # get parent gcf for bgc
            bgc_like_gcf = [gcf for gcf in self.gcfs if bgc_like.bgc_id in [b.bgc_id for b in gcf.bgcs]][0]
            # gather AS classes and convert to names in scoring dict
            as_classes = self.class_links.convert_as_classes(bgc_like.product_prediction.split("."))
            bgc_like_classes_dict = {"bigscape_class": bgc_like_gcf.bigscape_class,  # str - always one bigscape class right?
                                     "as_classes": as_classes}  # list(str)
        else:
            as_classes = self.class_links.convert_as_classes(self.class_links.get_gcf_as_classes(bgc_like, 0.5))
            bgc_like_classes_dict = {"bigscape_class": bgc_like.bigscape_class,  # str - always one bigscape class right?
                                     "as_classes": as_classes}  # list(str)
        if is_spectrum:
            # list of list of tuples/None - todo: add to spectrum object
            spec_like_classes = self.canopus.spectra_classes.get(str(spec_like.spectrum_id))
            spec_like_classes_names = self.canopus.spectra_classes_names
            spec_like_classes_names_inds = self.canopus.spectra_classes_names_inds
        else:  # molfam
            spec_like_classes = self.canopus.molfam_classes.get(str(spec_like.family_id))
            spec_like_classes_names = self.canopus.molfam_classes_names
            spec_like_classes_names_inds = self.canopus.molfam_classes_names_inds

        scores = []
        std_score = 0  # if link not recorded in scores (mibig) return this score
        # todo: calculate scores coming from spectrum side
        # loop through classes that are possible to link (names in link object)
        for bgc_class_name in self.class_links.bgc_class_names:
            if bgc_class_name == "mibig_classes":
                # treat specially as bigscape class needs to be translated to mibig class
                bigscape_class = bgc_like_classes_dict["bigscape_class"]
                # convert bigscape class to mibig class
                bgc_like_classes = [self.class_links.bigscape_mibig_conversion.get(bigscape_class)]
            else:
                bgc_like_classes = bgc_like_classes_dict.get(bgc_class_name)
            if bgc_like_classes and spec_like_classes:
                for bgc_class in bgc_like_classes:
                    for chem_class_name in self.class_links.chem_class_names:
                        # does info exist for this spectra class level
                        spec_class_i = spec_like_classes_names_inds.get(chem_class_name)
                        if spec_class_i:
                            spec_class_options = spec_like_classes[spec_class_i][0]  # get class - for now only use first
                            if spec_class_options:  # if there is a class at this lvl
                                spec_class = spec_class_options[0]  # is a tuple of (name, score) so take [0]
                                score = self.class_links.class_links[bgc_class_name][chem_class_name]\
                                    .get(bgc_class,{}).get(spec_class, std_score)
                                result_tuple = (score, bgc_class_name, chem_class_name, bgc_class, spec_class)
                                scores.append(result_tuple)
        return sorted(scores, reverse=True)

    def npclass_score(self, obj, target, method = 'main'):
        """Return sorted class link scores for scoring obj and target

        The input objects can be any mix of the following NPLinker types:
            - BGC
            - GCF
            - Spectrum
            - MolecularFamily

        Args:
            - obj: one of the possible input objects
            - target: one of the possible input objects
            - method: str, which classification method to use for spectra. default is 'main'
                options:
                    -'main': use the main method - currently canopus
                    -'mix': use main method first, when no main classes present, use the others
                    if present
                    -'canopus': use only canopus class predictions
                    -'molnetenhancer': use only molnetenhancer class predictions
        Returns:
            List of tuple of
            (score, obj_class_lvl, target_class_lvl, obj_class, target_class)
            (float, str, str, str, str)
            List will be empty if either BGC or spectrum classes are missing
        """
        # assess what is obj and target
        spec_like = obj
        bgc_like = target
        bgc_to_spec = False
        if isinstance(obj, BGC) or isinstance(obj, GCF):
            bgc_like = obj
            spec_like = target
            bgc_to_spec = True

        # assess if bgc or gcf and if spectrum or molfam
        is_spectrum = False
        is_bgc = False
        if isinstance(spec_like, Spectrum):
            is_spectrum = True
        if isinstance(bgc_like, BGC):
            is_bgc = True

        # assess method
        # todo: read options from NPLinker object - add option if the canopus/mne results are read correctly
        method_options = self.class_predict_options
        assert method in method_options, (
            f"NPClass method should be one of method options: {method_options}, if your method is not "+
            "in the options check if the class predictions (canopus, etc.) are loaded correctly")

        # gather correct classes based on input, dict for bgcs and list for spec
        if is_bgc:
            # get parent gcf for bgc
            bgc_like_gcf = [gcf for gcf in self.gcfs if bgc_like.bgc_id in [b.bgc_id for b in gcf.bgcs]][0]
            # gather AS classes and convert to names in scoring dict
            as_classes = self.class_links.convert_as_classes(bgc_like.product_prediction.split("."))
            bgc_like_classes_dict = {"bigscape_class": bgc_like_gcf.bigscape_class,  # str - always one bigscape class right?
                                     "as_classes": as_classes}  # list(str)
        else:
            as_classes = self.class_links.convert_as_classes(self.class_links.get_gcf_as_classes(bgc_like, 0.5))
            bgc_like_classes_dict = {"bigscape_class": bgc_like.bigscape_class,  # str - always one bigscape class right?
                                     "as_classes": as_classes}  # list(str)

        # gather classes for spectra, choose right method
        # choose the main method here by including it as 'main' in the method parameter
        use_canopus = 'main' in method or 'canopus' in method or 'mix' in method
        use_mne = 'molnetenhancer' in method or 'mix' in method
        spec_like_classes, spec_like_classes_names, spec_like_classes_names_inds = (None, None, None)
        # the order in which the classes are read, determines the priority (now: first canopus, then mne)
        if use_canopus and not spec_like_classes:
            if is_spectrum:
                # list of list of tuples/None - todo: add to spectrum object
                # take only 'best' (first) classification per ontology level
                all_classes = self.canopus.spectra_classes.get(str(spec_like.spectrum_id))
                if all_classes:
                    spec_like_classes = [cls_per_lvl for lvl in all_classes
                                         for i, cls_per_lvl in enumerate(lvl) if i==0]
                spec_like_classes_names = self.canopus.spectra_classes_names
                spec_like_classes_names_inds = self.canopus.spectra_classes_names_inds
            else:  # molfam
                all_classes = self.canopus.molfam_classes.get(str(spec_like.family_id))
                if all_classes:
                    spec_like_classes = [cls_per_lvl for lvl in all_classes
                                         for i, cls_per_lvl in enumerate(lvl) if i==0]
                spec_like_classes_names = self.canopus.molfam_classes_names
                spec_like_classes_names_inds = self.canopus.molfam_classes_names_inds
        if use_mne and not spec_like_classes:  # if mne or when main/canopus does not get classes
            if is_spectrum:
                spec_like_classes = self.molnetenhancer.spectra_classes(spec_like.spectrum_id)
            else: # molfam
                spec_like_classes = self.molnetenhancer.molfam_classes.get(str(spec_like.family_id))
            # classes are same for molfam and spectrum so names are irrespective of is_spectrum
            spec_like_classes_names = self.molnetenhancer.spectra_classes_names
            spec_like_classes_names_inds = self.molnetenhancer.spectra_classes_names_inds


        scores = []  # this will be returned if one of the class sides is absent
        std_score = 0  # if link not recorded in scores (mibig) return this score
        # loop through classes that are possible to link (names in class_link object)
        for bgc_class_name in self.class_links.bgc_class_names:
            if bgc_class_name == "mibig_classes":
                # treat specially as bigscape class needs to be translated to mibig class
                bigscape_class = bgc_like_classes_dict["bigscape_class"]
                # convert bigscape class to mibig class
                bgc_like_classes = [self.class_links.bigscape_mibig_conversion.get(bigscape_class)]
            else:
                bgc_like_classes = bgc_like_classes_dict.get(bgc_class_name)
            if bgc_like_classes and spec_like_classes:  # check for classes from both sides
                for bgc_class in bgc_like_classes:
                    for chem_class_name in self.class_links.chem_class_names:
                        # does info exist for this spectrum class level, return index for class level
                        spec_class_level_i = spec_like_classes_names_inds.get(chem_class_name)
                        if spec_class_level_i:
                            spec_class_tup = spec_like_classes[spec_class_level_i]
                            if spec_class_tup:  # if there is a class at this lvl
                                spec_class = spec_class_tup[0]  # is a tuple of (name, score) so take [0]

                                # determine direction of scoring: BGC -> spectrum
                                if bgc_to_spec:
                                    score = self.class_links.class_links[bgc_class_name][chem_class_name]\
                                        .get(bgc_class,{}).get(spec_class, std_score)
                                    result_tuple = (score, bgc_class_name, chem_class_name, bgc_class, spec_class)
                                else:  # spectrum -> BGC
                                    score = self.class_links.class_links[chem_class_name][bgc_class_name]\
                                        .get(spec_class, {}).get(bgc_class, std_score)
                                    result_tuple = (score, chem_class_name, bgc_class_name, spec_class, bgc_class)
                                scores.append(result_tuple)
        return sorted(scores, reverse=True)

    @property
    def class_links(self):
        return self._class_links

    @property
    def canopus(self):
        return self._canopus

    @property
    def molnetenhancer(self):
        return self._molnetenhancer

    @property
    def class_predict_options(self):
        return self._class_predict_options
