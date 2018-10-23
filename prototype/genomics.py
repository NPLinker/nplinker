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


def loadBGC_from_cluster_files(network_file_list,ann_file_list):
    strain_id_dict = {}
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
                try:
                    try:
                        strain_name = strain_id_dict[name.split('_')[0]]
                    except:
                        strain_name = strain_id_dict[name.split('.')[0]]
                except:
                    # it's a MiBIG one
                    strain_name = 'MiBIG'
                print name,strain_name


    return 1,2,3

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
    

