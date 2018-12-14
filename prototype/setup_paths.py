import os,glob

def setup_paths(root_dir):
    MIBIG_JSON_DIR = os.path.join(root_dir,'mibig_json')
    print("MIBIG_JSON_DIR:",MIBIG_JSON_DIR)

    spec_dir = os.path.join(root_dir,'spectra')
    MGF_FILE = glob.glob(os.path.join(spec_dir,'*.mgf'))[0]

    print("MGF_FILE:",MGF_FILE)
    ANNOTATION_FILES = glob.glob(os.path.join(spec_dir,'*.annotations.*'))
    print("ANNOTATION_FILES:"," ".join(ANNOTATION_FILES))

    EDGES_FILE = glob.glob(os.path.join(spec_dir,"*.pairsinfo"))[0]
    print("EDGES_FILE:",EDGES_FILE)

    NODES_FILE = glob.glob(os.path.join(spec_dir,'reformated.tsv'))
    print("NODES_FILE:",NODES_FILE)

    ROOT_PATH = os.path.join(root_dir,'bigscape')
    ANTISMASH_DIR = os.path.join(root_dir,'antismash')

    return MIBIG_JSON_DIR,NODES_FILE,ANNOTATION_FILES,MGF_FILE,EDGES_FILE,ROOT_PATH,ANTISMASH_DIR


def setup_paths_simon():




    MIBIG_JSON_DIR = "/Users/simon/Desktop/carnegie_bigscape/mibig/mibig_json-1.4/"
    # PATH_MS2LDA = DATASET

    # annotations
    NODES_FILE = "/Users/simon/Dropbox/BioResearch/Meta_clustering/MS2LDA/carnegie_gnps_data/ProteoSAFe-METABOLOMICS-SNETS-2cdc2aa3-download_clustered_spectra-2/clusterinfosummarygroup_attributes_withIDs_withcomponentID/reformated_derep.tsv"
    DEREP_ANNOTATIONS = "/Users/simon/Dropbox/BioResearch/Meta_clustering/MS2LDA/carnegie_gnps_data/cluster_to_derep.tsv"

    MGF_FILE = "/Users/simon/Dropbox/BioResearch/Meta_clustering/MS2LDA/carnegie_gnps_data/ProteoSAFe-METABOLOMICS-SNETS-2cdc2aa3-download_clustered_spectra-2/METABOLOMICS-SNETS-2cdc2aa3-download_clustered_spectra-main.mgf"
    EDGES_FILE = "/Users/simon/Dropbox/BioResearch/Meta_clustering/MS2LDA/carnegie_gnps_data/ProteoSAFe-METABOLOMICS-SNETS-2cdc2aa3-download_clustered_spectra-2/networkedges_selfloop/1eca77da22c84ae1aee490d0e4bc26f3.pairsinfo"

    FOLDERS = ['NRPS','Others','PKSI','PKS-NRP_Hybrids','PKSother','RiPPs','Saccharides','Terpene']
    FOLDERS = ['NRPS']
    # GCF file path
    ROOT_PATH = "/Users/simon/Desktop/carnegie_bigscape/bigscape-carnegie-glocal/network_files/2018-11-16_12-01-18_hybrids_glocal/"
    ANTISMASH_DIR = "/Users/simon/Desktop/carnegie_bigscape/antismash/"
    return MIBIG_JSON_DIR,NODES_FILE,DEREP_ANNOTATIONS,MGF_FILE,EDGES_FILE,FOLDERS,ROOT_PATH,ANTISMASH_DIR