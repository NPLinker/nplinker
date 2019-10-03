import csv, glob, os, json

import numpy as np

import aa_pred
from genomics_utilities import get_known_cluster_blast
from genomics_utilities import get_smiles

from logconfig import LogConfig
logger = LogConfig.getLogger(__file__)

class BGC(object):
    def __init__(self, id, strain, name, bigscape_class, product_prediction, description=None):
        self.id = id
        self.strain = strain
        self.name = name
        self.bigscape_class = bigscape_class
        self.product_prediction = product_prediction
        self.parent = None
        self.description = description

        self.antismash_file = None
        self._aa_predictions = None
        self._known_cluster_blast = None
        self._smiles = None
        self._smiles_parsed = False

        self.edges = []

    def __repr__(self):
        return str(self)

    def __str__(self):
        return self.__class__.__name__ + "(name=" + self.name + ", strain=" + str(self.strain) + ")"
    
    @property
    def smiles(self):
        if self._smiles is not None or self._smiles_parsed: 
            return self._smiles

        if self.antismash_file is None:
            return None 
    
        self._smiles = get_smiles(self)
        self._smiles_parsed = True
        logger.debug('SMILES for {} = {}'.format(self, self._smiles))
        return self._smiles

    @property
    def aa_predictions(self):
        # Load aa predictions and cache them
        self._aa_predictions = None
        if self._aa_predictions is None:
            self._aa_predictions = {}
            if self.antismash_file is not None:
                for p in aa_pred.predict_aa(self.antismash_file):
                    self._aa_predictions[p[0]] = p[1]
        return [self._aa_predictions]

    @property
    def known_cluster_blast(self):
        if self._known_cluster_blast is None:
            self._known_cluster_blast = get_known_cluster_blast(self)
        return self._known_cluster_blast


class GCF(object):
    def __init__(self, gcf_id, gnps_class):
        self.id = -1
        self.gcf_id = gcf_id
        self.gnps_class = gnps_class
        try:
            self.short_gcf_id = gcf_id.split(os.sep)[-1]
        except:
            self.short_gcf_id = self.gcf_id
        self.bgc_list = []
        self.random_gcf = None

        self._aa_predictions = None
        self.strains_lookup = None

    def __str__(self):
        return 'GCF(id={}, class={}, gcf_id={})'.format(self.id, self.gnps_class, self.short_gcf_id)

    def __repr__(self):
        return str(self)

    def add_bgc(self, bgc):
        self.bgc_list.append(bgc)
        bgc.parent = self

    @property
    def strains(self):
        return [bgc.strain for bgc in self.bgc_list]

    def only_mibig(self):
        return len(self.bgc_list) == self.num_mibig_bgcs

    def has_mibig(self):
        return self.num_mibig_bgcs > 0

    def has_strain(self, strain):
        for bgc in self.bgc_list:
            if bgc.strain == strain:
                return True
        return False

    def add_random(self, bgc_list):
        self.random_gcf = RandomGCF(self, bgc_list)

    @property
    def non_mibig_bgcs(self):
        return list(filter(lambda bgc: not isinstance(bgc, MiBIGBGC)))

    @property
    def mibig_bgcs(self):
        return list(filter(lambda bgc: isinstance(bgc, MiBIGBGC)))

    @property
    def num_mibig_bgcs(self):
        return len(self.mibig_bgcs)

    @property
    def num_non_mibig_bgcs(self):
        return len(self.bgc_list) - self.num_mibig_bgcs

    @property
    def aa_predictions(self):
        """
        Return the predicted AAs for the GCF
        """
        if self._aa_predictions is None:
            bgc_aa_prob = []
            for bgc_count, bgc in enumerate(self.bgc_list):
                bgc_aa_prob.extend(bgc.aa_predictions)
            self._aa_predictions = bgc_aa_prob

        return self._aa_predictions


class RandomGCF(object):

    def __init__(self, real_gcf, bgc_list):
        self.real_gcf = real_gcf
        n_bgc = self.real_gcf.num_non_mibig_bgcs
        print(n_bgc)
        # select n_bgc bgcs from the bgc_list
        # (convert back to list for consistency with GCF)
        self.bgc_list = list(np.random.choice(bgc_list, n_bgc, replace=False))

    @property
    def strains(self):
        return [bgc.strain for bgc in self.bgc_list]

    def has_strain(self, strain):
        for bgc in self.bgc_list:
            if bgc.strain == strain:
                return True
        return False

class MiBIGBGC(BGC):

    def __init__(self, id, name, product_prediction):
        super(MiBIGBGC, self).__init__(id, name, name, None, product_prediction)

def loadBGC_from_cluster_files(cluster_file_list, ann_file_list, network_file_list, antismash_dir, antismash_filenames, antismash_format='default', mibig_bgc_dict=None):
    strain_id_dict = {}
    strain_dict = {}
    gcf_dict = {}
    gcf_list = []
    strain_list = []
    bgc_list = []

    # TODO might need to change this later depending on packaging etc
    with open(os.path.join(os.path.dirname(__file__), 'strain_ids.csv'), 'r') as f:
        reader = csv.reader(f)
        for line in reader:
            strain_id_dict[line[0]] = line[1]

    metadata = {}
    # parse the annotation files (<dataset>/bigscape/<cluster_name>/Network_Annotations_<cluster_name>.tsv
    # these contain fields:
    # - BGC name/ID [0]
    # - "Accession ID" [1]
    # - Description [2]
    # - Product prediction [3]
    # - Bigscape class [4]
    # - Organism [5]
    # - Taxonomy [6]
    for a in ann_file_list:
        with open(a, 'rU') as f:
            reader = csv.reader(f, delimiter='\t')
            heads = next(reader)
            for line in reader:
                metadata[line[0]] = line

    num_mibig = 0
    num_missing_antismash = 0

    bgc_lookup = {}

    # "cluster files" are the various <class>_clustering_c0.xx.tsv files
    # - BGC name
    # - cluster ID
    for filename in cluster_file_list:
        gnps_class = os.path.split(filename)[-1]
        gnps_class = gnps_class[:gnps_class.index('_')]
        with open(filename, 'rU') as f:
            reader = csv.reader(f, delimiter='\t')
            heads = next(reader)
            for line in reader:
                name = line[0]
                #family = os.path.split(filename)[-1] + ":" + line[1]
                family_id = int(line[1])
                if name.startswith("BGC"):
                    strain_name = 'MiBIG'
                else:
                    # TODO is this doing the correct thing in all cases??
                    try:
                        try:
                            strain_name = strain_id_dict[name.split('_')[0]]
                        except:
                            strain_name = strain_id_dict[name.split('.')[0]]
                    except:
                        # logger.warning("strain lookup failed for '%s'" % name)
                        if '_' in name:
                            strain_name = name.split('_')[0] 
                        else:
                            strain_name = name.split('.')[0]

                if strain_name not in strain_dict:
                    new_strain = strain_name
                    strain_dict[strain_name] = new_strain
                    strain_list.append(new_strain)

                strain = strain_dict[strain_name]
                tokens = name.split('.')
                clusterid = tokens[-1]
                rest = '.'.join(tokens[:-1])

                metadata_line = metadata[name]
                description = metadata_line[2]
                bigscape_class = metadata_line[4]
                product_prediction = metadata_line[3]

                # make a BGC object, reusing existing objects if they represent the same physical thing
                if not strain_name == 'MiBIG':
                    new_bgc = bgc_lookup.get(name, BGC(len(bgc_list), strain, name, bigscape_class, product_prediction, description))
                    if antismash_dir:
                        if antismash_format == 'flat':
                            antismash_filename = os.path.join(antismash_dir, new_bgc.name + '.gbk')
                            if not os.path.exists(antismash_filename):
                                logger.warn('!!! Missing antismash file: {}'.format(antismash_filename))
                                num_missing_antismash += 1
                                # return None, None, None
                            new_bgc.antismash_file = antismash_filename
                        else:
                            new_bgc.antismash_file = antismash_filenames.get(new_bgc.name, None)
                            if new_bgc.antismash_file is None:
                                logger.warning('Failed to find an antiSMASH file for {}'.format(new_bgc.name))
                                num_missing_antismash += 1
                    bgc_list.append(new_bgc)
                else:
                    num_mibig += 1
                    # TODO should this not attempt to create any MiBIGBGC's if the dict
                    # is not supplied??
                    # TODO any reason not to supply the metadata fields that aren't set by
                    # make_mibig_bgc_dict since metadata_line is available here?
                    if mibig_bgc_dict is not None:
                        try:
                            new_bgc = mibig_bgc_dict[name.split('.')[0]]
                        except KeyError:
                            new_bgc = MiBIGBGC(len(bgc_list), name, product_prediction)

                    new_bgc.bigscape_class = bigscape_class
                    new_bgc.description = description
                    # TODO should add to bgc_list too???
                    bgc_list.append(new_bgc)

                if family_id not in gcf_dict:
                    new_gcf = GCF(family_id, gnps_class)
                    gcf_dict[family_id] = new_gcf
                    gcf_list.append(new_gcf)
                gcf_dict[family_id].add_bgc(new_bgc)

                bgc_lookup[new_bgc.name] = new_bgc
    
    if num_missing_antismash > 0:
        logger.warn('{}/{} antiSMASH files could not be found!'.format(num_missing_antismash, len(bgc_list)))
        # print(list(antismash_filenames.keys()))

    logger.debug('Loading .network files')
    for filename in network_file_list:
        with open(filename, 'rU') as f:
            reader = csv.reader(f, delimiter='\t')
            heads = next(reader)
            # try to look up bgc IDs
            for line in reader:
                bgc_src = bgc_lookup[line[0]]
                bgc_dst = bgc_lookup[line[1]]

                bgc_src.edges.append(bgc_dst.id)

    print('# MiBIG BGCs = {}, non-MiBIG BGCS = {}, total bgcs = {}, GCFs = {}'.format(
                        num_mibig, len(bgc_list) - num_mibig, len(bgc_list), len(gcf_dict)))
    # Assign unique ids (int)
    # TODO do this above
    for i, gcf in enumerate(gcf_list):
        gcf.id = i
    return gcf_list, bgc_list, strain_list

# this is really slow. But hey, it works. Finally.
def find_antismash_file_flat(antismash_dir, bgc_name):
    all_gbk_files = glob.glob(antismash_dir + os.sep + '*.gbk')
    last_bit = [o.split(os.sep)[-1] for o in all_gbk_files]
    bgc_file_name = bgc_name + '.gbk'
    if bgc_file_name in last_bit:
        idx = last_bit.index(bgc_file_name)
        return all_gbk_files[idx]
    else:
        print("NOOO",bgc_name)
        return None

def find_antismash_file_old(antismash_dir, bgc_name):
    subdirs = [s.split(os.sep)[-1] for s in glob.glob(antismash_dir + os.sep+'*')]
    if bgc_name.startswith('BGC'):
        print("No file for MiBIG BGC")
        return None # MiBIG BGC

    # this code is nasty... :-)
    name_tokens = bgc_name.split('_')
    found = False
    for i in range(len(name_tokens)):
        sub_name = '_'.join(name_tokens[:i])
        if sub_name in subdirs:
            found = True
            found_name = sub_name
    if not found:
        name_tokens = bgc_name.split('.')[0]
        for i in range(len(name_tokens)):
            sub_name = '.'.join(name_tokens[:i])
            if sub_name in subdirs:
                found = True
                found_name = sub_name
    if not found:
        print("Can't find antiSMASH info for ", bgc_name)
        return None
#     print found_name
    dir_contents = glob.glob(antismash_dir + os.sep + found_name + os.sep + '*.gbk')
    cluster_names = [d.split('.')[-2] for d in dir_contents]
    this_name = bgc_name.split('.')[-1]
    try:
        antismash_name = dir_contents[cluster_names.index(bgc_name.split('.')[-1])]
    except:
        print(bgc_name)
        print(cluster_names)
        print()
        print()
        return None
    return antismash_name

#
# not currently used!
#
def loadBGC_from_node_files(file_list):
    strain_id_dict = {}
    with open('strain_ids.csv', 'r') as f:
        reader = csv.reader(f)
        for line in reader:
            strain_id_dict[line[0]] = line[1]
    
    strain_dict = {}
    gcf_dict = {}
    bgc_list = []
    gcf_list = []
    strain_list = []
    for filename in file_list:
        with open(filename, 'rU') as f:
            reader = csv.reader(f)
            heads = next(reader)
            name_pos = heads.index("shared name")
            description_pos = heads.index("Description")
            bigscape_class_pos = heads.index("BiG-SCAPE class")
            product_prediction_pos = heads.index("Product Prediction")
            family_pos = heads.index("Family Number")
            for line in reader:
                name = line[name_pos]
                try:
                    try:
                        strain_name = strain_id_dict[name.split('_')[0]]
                    except:
                        strain_name = strain_id_dict[name.split('.')[0]]
                except:
                    # it's a MiBIG one
                    strain_name = 'MiBIG'
                if strain_name not in strain_dict:
                    new_strain = strain_name
                    strain_dict[strain_name] = new_strain
                    strain_list.append(new_strain)
                strain = strain_dict[strain_name]


                tokens = name.split('.')
                clusterid = tokens[-1]
                rest = '.'.join(tokens[:-1])
                # print name,rest,clusterid
                description = line[description_pos]
                bigscape_class = line[bigscape_class_pos]
                product_prediction = line[product_prediction_pos]
                family = filename + " " + line[family_pos]

                # make a BGC object
                if not strain_name == 'MiBIG':
                    new_bgc = BGC(strain, name, bigscape_class, product_prediction)
                else:
                    new_bgc = MiBIGBGC(name, product_prediction)
                bgc_list.append(new_bgc)

                if family not in gcf_dict:
                    new_gcf = GCF(family)
                    gcf_dict[family] = new_gcf
                    gcf_list.append(new_gcf)
                gcf_dict[family].add_bgc(new_bgc)

    return gcf_list, bgc_list, strain_list

def load_mibig_map(filename='mibig_gnps_links_q3_loose.csv'):
    mibig_map = {}
    with open(filename, 'rU') as f:
        reader = csv.reader(f)
        heads = next(reader)

        for line in reader:
            bgc = line[0]
            if bgc in mibig_map:
                mibig_map[bgc].append(line[3])
            else:
                mibig_map[bgc] = [line[3]]
    return mibig_map


def load_mibig_library_json(mibig_json_directory):
    mibig = {}
    files = glob.glob(mibig_json_directory + os.sep + '*.json')
    print("Found {} MiBIG json files".format(len(files)))
    for file in files:
        with open(file, 'r') as f:
            bgc_id = file.split(os.sep)[-1].split('.')[0]
            mibig[bgc_id] = json.load(f)
    return mibig

# TODO this will need to have consistent bgc.id values if called
# separately during the loading process, and/or have existing
# MiBIGBGC objects merged with this set...
def make_mibig_bgc_dict(mibig_json_directory):
    mibig_dict = load_mibig_library_json(mibig_json_directory)
    mibig_bgc_dict = {}
    i = 0
    for name, data in list(mibig_dict.items()):
        new_bgc = MiBIGBGC(i, data['general_params']['mibig_accession'], data['general_params']['biosyn_class'])
        mibig_bgc_dict[data['general_params']['mibig_accession']] = new_bgc
        i += 1
    return mibig_bgc_dict
