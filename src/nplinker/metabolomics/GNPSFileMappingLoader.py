import csv
from nplinker.metabolomics.IFileMappingLoader import IFileMappingLoader
from nplinker.metabolomics.load_gnps import identify_gnps_format


class GNPSFileMappingLoader(IFileMappingLoader):

    def __init__(self, filename):
        self._mapping = load_file_mappings(filename)

def load_file_mappings(filename):
    with open(filename, mode='rt', encoding='utf-8') as file:
        csv.reader(file, delimiter='\t')