import csv, glob, os, json

import numpy as np

from .aa_pred import predict_aa
from .genomics_utilities import get_smiles

from .strains import Strain, StrainCollection

from .logconfig import LogConfig
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

        self.edges = set()

    def __repr__(self):
        return str(self)

    def __str__(self):
        return '{}(id={}, name={}, strain={})'.format(self.__class__.__name__, self.id, self.name, self.strain)

    def __eq__(self, other):
        return self.id == other.id

    def __hash__(self):
        return self.id

    @property
    def is_mibig(self):
        return self.name.startswith('BGC')

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
                for p in predict_aa(self.antismash_file):
                    self._aa_predictions[p[0]] = p[1]
        return [self._aa_predictions]

class MiBIGBGC(BGC):

    def __init__(self, id, strain, name, product_prediction):
        super(MiBIGBGC, self).__init__(id, strain, name, None, product_prediction)

    def __eq__(self, other):
        return self.id == other.id

    def __hash__(self):
        return self.id


class GCF(object):
    def __init__(self, id, gcf_id, product_type):
        self.id = id
        self.gcf_id = gcf_id
        self.product_type = product_type
        self.bgcs = set()
        self.classes = set()

        self._aa_predictions = None
        self.strains = StrainCollection()
        self.strains_lookup = {}
        self.dataset_strains = None 

    def __str__(self):
        return 'GCF(id={}, class={}, gcf_id={}, strains={})'.format(self.id, self.product_type, self.gcf_id, len(self.strains))

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return self.id == other.id

    def __hash__(self):
        return self.id

    def add_bgc(self, bgc):
        self.bgcs.add(bgc)
        self.classes.add(bgc.bigscape_class)
        # TODO possible to have multiple parents??
        bgc.parent = self
        self.strains.add(bgc.strain)
        self.strains_lookup[bgc.strain] = bgc

    def has_strain(self, strain):
        return strain in self.strains

    def bgc_for_strain(self, strain):
        return self.strains_lookup[strain]

    def only_mibig(self):
        return len(self.bgcs) == self.num_mibig_bgcs

    def has_mibig(self):
        return self.num_mibig_bgcs > 0

    @property
    def non_mibig_bgcs(self):
        return list(filter(lambda bgc: not isinstance(bgc, MiBIGBGC), self.bgcs))

    @property
    def mibig_bgcs(self):
        return list(filter(lambda bgc: isinstance(bgc, MiBIGBGC), self.bgcs))

    @property
    def num_mibig_bgcs(self):
        return len(self.mibig_bgcs)

    @property
    def num_non_mibig_bgcs(self):
        return len(self.bgcs) - self.num_mibig_bgcs

    @property
    def aa_predictions(self):
        """
        Return the predicted AAs for the GCF
        """
        if self._aa_predictions is None:
            bgc_aa_prob = []
            for bgc_count, bgc in enumerate(self.bgcs):
                if not bgc.name.startswith('BGC'):
                    bgc_aa_prob.extend(bgc.aa_predictions)
            self._aa_predictions = bgc_aa_prob

        return self._aa_predictions

def loadBGC_from_cluster_files(strains, cluster_file_dict, ann_file_dict, network_file_dict, mibig_bgc_dict, antismash_dir, antismash_filenames, antismash_format, antismash_delimiters):
    gcf_dict = {}
    gcf_list = []
    metadata = {}

    # parse the annotation files (<dataset>/bigscape/<cluster_name>/Network_Annotations_<cluster_name>.tsv
    # these contain fields:
    # - BGC name/ID [0]
    # - "Accession ID" [1]
    # - Description [2]
    # - Product prediction [3]
    # - Bigscape product type/class [4]
    # - Organism [5]
    # - Taxonomy [6]
    for a in ann_file_dict.values():
        with open(a, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader) # skip headers
            for line in reader:
                metadata[line[0]] = line

    num_mibig = 0
    num_missing_antismash = 0

    bgc_lookup = {}
    internal_bgc_id = len(mibig_bgc_dict) # start numbering BGCs from here
    internal_gcf_id = 0

    bgc_list = [v for v in mibig_bgc_dict.values()]

    unknown_strains = {}

    logger.info('Using antiSMASH filename delimiters {}'.format(antismash_delimiters))

    # "cluster files" are the various <class>_clustering_c0.xx.tsv files
    # - BGC name
    # - cluster ID
    for product_type, filename in cluster_file_dict.items():
        product_type = os.path.split(filename)[-1]
        product_type = product_type[:product_type.index('_')]
        with open(filename, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader) # skip headers
            for line in reader:
                name = line[0]
                family_id = int(line[1])
                if name.startswith('BGC'):
                    # removing the .<digit> suffix
                    nname = name[:name.index('.')]
                    strain = strains.lookup(nname)
                    if strain is None:
                        raise Exception('Unknown MiBIG BGC: original={} / parsed={}'.format(name, nname))
                else:
                    parsednames = [name[:name.index(d)] for d in antismash_delimiters if name.find(d) != -1]
                    found = False
                    for parsedname in parsednames:
                        strain = strains.lookup(parsedname)
                        if strain is not None:
                            found = True
                            break
                        else:
                            # TODO hack to get crusemann working, should really update strain mappings?
                            for i in range(1, 3, 1):
                                tmp = '{}.{}'.format(parsedname, i)
                                strain = strains.lookup(tmp)
                                if strain is not None:
                                    found = True
                                    break
                            if found:
                                break

                    if not found:
                        logger.warning('Unknown strain ID: {} (from file {})'.format(name, filename))
                        unknown_strains[name] = filename
                        continue

                metadata_line = metadata[name]
                description = metadata_line[2]
                bigscape_class = metadata_line[4]
                product_prediction = metadata_line[3]

                # make a BGC object, reusing existing objects if they represent the same physical thing
                if not strain.id.startswith('BGC'):
                    if name in bgc_lookup:
                        new_bgc = bgc_lookup.get(name)
                    else:
                        # create a new BGC, increment internal ID and add to the list
                        new_bgc = BGC(internal_bgc_id, strain, name, bigscape_class, product_prediction, description)
                        internal_bgc_id += 1
                        bgc_list.append(new_bgc)

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
                            # TODO should this actually search for all possible names based on strain
                            # aliases instead of just this name? 
                            if new_bgc.antismash_file is None:
                                logger.warning('Failed to find an antiSMASH file for {}'.format(new_bgc.name))
                                num_missing_antismash += 1
                else:
                    num_mibig += 1
                    # TODO any reason not to supply the metadata fields that aren't set by
                    # make_mibig_bgc_dict since metadata_line is available here?
                    try:
                        new_bgc = mibig_bgc_dict[strain.id]
                    except KeyError:
                        raise Exception('Unknown MiBIG: {}'.format(strain.id))

                    new_bgc.bigscape_class = bigscape_class
                    new_bgc.description = description

                if family_id not in gcf_dict:
                    new_gcf = GCF(internal_gcf_id, family_id, bigscape_class)
                    gcf_dict[family_id] = new_gcf
                    gcf_list.append(new_gcf)
                    internal_gcf_id += 1

                gcf_dict[family_id].add_bgc(new_bgc)

                bgc_lookup[new_bgc.name] = new_bgc
    
    if num_missing_antismash > 0:
        logger.warn('{}/{} antiSMASH files could not be found!'.format(num_missing_antismash, len(bgc_list)))
        # print(list(antismash_filenames.keys()))

    logger.info('# MiBIG BGCs = {}, non-MiBIG BGCS = {}, total bgcs = {}, GCFs = {}, strains={}'.format(
                        num_mibig, len(bgc_list) - num_mibig, len(bgc_list), len(gcf_dict), len(strains)))

    # filter out irrelevant MiBIG BGCs (and MiBIG-only GCFs)
    bgc_list, gcf_list, strains = filter_mibig_bgcs(bgc_list, gcf_list, strains)
    # update lookup table as well
    bgc_lookup = {bgc.name: bgc for bgc in bgc_list}

    logger.info('# after filtering, total bgcs = {}, GCFs = {}, strains={}'.format(len(bgc_list), len(gcf_list), len(strains)))

    # load edge info - note that this should be done AFTER the filtering step above
    # so that it won't leave us with edges for BGCs that are no longer present
    logger.debug('Loading .network files')
    for filename in network_file_dict.values():
        with open(filename, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader) # skip headers
            # try to look up bgc IDs
            for line in reader:
                for i in range(2):
                    if line[i].startswith('BGC'):
                        # removing the .<digit> suffix
                        line[i] = line[i][:line[i].index('.')]

                if line[0] not in bgc_lookup or line[1] not in bgc_lookup:
                    # should indicate that one or both of these BGCs have been filtered out above
                    continue

                bgc_src = bgc_lookup[line[0]]
                bgc_dst = bgc_lookup[line[1]]
                bgc_src.edges.add(bgc_dst.id)

    return gcf_list, bgc_list, strains, unknown_strains

def filter_mibig_bgcs(bgcs, gcfs, strains):
    # remove the following MiBIG BGCs:
    # - parent attr is None (indicating never added to a GCF)
    # - any instances in a GCF with no other non-MiBIG BGCs
    to_remove_gcfs = set()
    to_remove_bgcs = set()
    for gcf in gcfs:
        # can ignore any with no MiBIGs
        if gcf.num_mibig_bgcs == 0:
            continue
        # now know this GCF has >0 MiBIG BGCs. Next step is to
        # check number of non-MiBIG BGCs, and if this is not at 
        # least 1, throw away both GCF and BGC(s)
        if gcf.num_non_mibig_bgcs == 0:
            to_remove_gcfs.add(gcf)
            for bgc in gcf.bgcs:
                to_remove_bgcs.add(bgc)
                strains.remove(bgc.strain)
                bgc.parent = None

    for bgc in bgcs:
        if bgc.parent is None:
            strains.remove(bgc.strain)

    logger.info('Filtering MiBIG BGCs: removing {} GCFs and {} BGCs'.format(len(to_remove_gcfs), len(to_remove_bgcs)))

    # for GCFs just remove those that appear in to_remove_gcfs
    new_gcf_list = [gcf for gcf in gcfs if gcf not in to_remove_gcfs]
    # for BGCs do similar but also get rid of the objects never added to a GCF in the first place
    new_bgc_list = [bgc for bgc in bgcs if bgc not in to_remove_bgcs and bgc.parent is not None]

    # keep internal IDs consecutive 
    for i in range(len(new_bgc_list)):
        new_bgc_list[i].id = i
        if new_bgc_list[i].parent is None:
            raise Exception(new_bgc_list[i])

    for i in range(len(new_gcf_list)):
        new_gcf_list[i].id = i

    return new_bgc_list, new_gcf_list, strains

def load_mibig_map(filename='mibig_gnps_links_q3_loose.csv'):
    mibig_map = {}
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        next(reader) # skip headers

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
    logger.info("Found {} MiBIG json files".format(len(files)))
    for file in files:
        with open(file, 'r') as f:
            bgc_id = file.split(os.sep)[-1].split('.')[0]
            mibig[bgc_id] = json.load(f)
    return mibig

def make_mibig_bgc_dict(strains, mibig_json_directory):
    mibig_dict = load_mibig_library_json(mibig_json_directory)
    mibig_bgc_dict = {}
    i = 0
    for name, data in list(mibig_dict.items()):
        accession = data['general_params']['mibig_accession']
        biosyn_class = data['general_params']['biosyn_class'][0]
        strain = Strain(accession)
        new_bgc = MiBIGBGC(i, strain, accession, biosyn_class)
        mibig_bgc_dict[accession] = new_bgc
        strains.add(strain)
        i += 1
    return mibig_bgc_dict
