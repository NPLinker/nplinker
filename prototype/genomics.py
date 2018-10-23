import csv,glob
import numpy as np

class Strain(object):
    def __init__(self,name):
        self.name = name
    def __str__(self):
        return self.name

class BGC(object):
    def __init__(self,strain,name,bigscape_class,product_prediction):
        self.strain = strain
        self.name = name
        self.bigscape_class = bigscape_class
        self.product_prediction = product_prediction
    def __str__(self):
        return self.name + "(" + str(self.strain) + ")"

class GCF(object):
    def __init__(self,gcf_id):
        self.gcf_id = gcf_id
        self.bgc_list = []
        self.random_gcf = None

    def add_bgc(self,bgc):
        self.bgc_list.append(bgc)

    def has_strain(self,strain):
        for bgc in self.bgc_list:
            if bgc.strain == strain:
                return True
        return False

    def add_random(self,strain_list):
        self.random_gcf = RandomGCF(self,strain_list)

class RandomGCF(object):
    def __init__(self,real_gcf,strain_list):
        n_strains = 0
        self.real_gcf = real_gcf
        for s in strain_list:
            if self.real_gcf.has_strain(s):
                n_strains += 1
        # select n_strains from strain_list
        self.strain_set = set(np.random.choice(strain_list,n_strains,replace = False))
    def has_strain(self,strain):
        if strain in self.strain_set:
            return True
        else:
            return False


class MiBIGBGC(BGC):
    def __init__(self,name,product_prediction):
        super(MiBIGBGC,self).__init__(None,name,None,product_prediction)


def loadBGC_from_cluster_files(network_file_list,ann_file_list,antismash_dir = None):
    strain_id_dict = {}
    strain_dict = {}
    gcf_dict = {}
    gcf_list = []
    strain_list = []
    bgc_list = []
    with open('strain_ids.csv','r') as f:
        reader = csv.reader(f)
        for line in reader:
            strain_id_dict[line[0]] = line[1]
    metadata = {}
    for a in ann_file_list:
        with open(a,'rU') as f:
            reader =  csv.reader(f,delimiter = '\t')
            heads = reader.next()
            for line in reader:
                metadata[line[0]] = line

    for filename in network_file_list:
        with open(filename,'rU') as f:
            reader = csv.reader(f,delimiter = '\t')
            heads = reader.next()
            for line in reader:
                name = line[0]
                family = filename + ":" + line[1]
                if name.startswith("BGC"):
                    strain_name = 'MiBIG'
                else:
                    try:
                        try:
                            strain_name = strain_id_dict[name.split('_')[0]]
                        except:
                            strain_name = strain_id_dict[name.split('.')[0]]
                    except:
                        print "NO STRAIN"
                if not strain_name in strain_dict:
                    new_strain = Strain(strain_name)
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

                # make a BGC object
                # the same BGC objects might be made more than once
                # because they appear in multiple clusterings
                if not strain_name == 'MiBIG':
                    new_bgc = BGC(strain,name,bigscape_class,product_prediction)
                    if antismash_dir:
                        new_bgc.antismash_file = find_antismash_file(antismash_dir,new_bgc.name)
                else:
                    new_bgc = MiBIGBGC(name,product_prediction)
                bgc_list.append(new_bgc)

                if not family in gcf_dict:
                    new_gcf = GCF(family)
                    gcf_dict[family] = new_gcf
                    gcf_list.append(new_gcf)
                gcf_dict[family].add_bgc(new_bgc)

    return gcf_list,bgc_list,strain_list


def find_antismash_file(antismash_dir,bgc_name):
    import glob,os
    subdirs = [s.split(os.sep)[-1] for s in glob.glob(antismash_dir + os.sep+'*')]
    if bgc_name.startswith('BGC'):
        print "No file for MiBIG BGC"
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
        print "Can't find antiSMASH info for ",bgc_name
        return None
#     print found_name
    dir_contents = glob.glob(antismash_dir + os.sep + found_name + os.sep + '*.gbk')
    cluster_names = [d.split('.')[-2] for d in dir_contents]
    this_name = bgc_name.split('.')[-1]
    try:
        antismash_name = dir_contents[cluster_names.index(bgc_name.split('.')[-1])]
    except:
        print bgc_name
        print cluster_names
        print
        print
    return antismash_name

def loadBGC_from_node_files(file_list):
    strain_id_dict = {}
    with open('strain_ids.csv','r') as f:
        reader = csv.reader(f)
        for line in reader:
            strain_id_dict[line[0]] = line[1]
    
    strain_dict = {}
    gcf_dict = {}
    bgc_list = []
    gcf_list = []
    strain_list = []
    for filename in file_list:
        with open(filename,'rU') as f:
            reader = csv.reader(f)
            heads = reader.next()
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
                if not strain_name in strain_dict:
                    new_strain = Strain(strain_name)
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
                    new_bgc = BGC(strain,name,bigscape_class,product_prediction)
                else:
                    new_bgc = MiBIGBGC(name,product_prediction)
                bgc_list.append(new_bgc)

                if not family in gcf_dict:
                    new_gcf = GCF(family)
                    gcf_dict[family] = new_gcf
                    gcf_list.append(new_gcf)
                gcf_dict[family].add_bgc(new_bgc)


                
    return gcf_list,bgc_list,strain_list
    

