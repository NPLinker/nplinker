import csv
from typing import Dict, Literal
from nplinker.metabolomics.IFileMappingLoader import IFileMappingLoader
from nplinker.metabolomics.load_gnps import identify_gnps_format


class GNPSFileMappingLoader(IFileMappingLoader):

    def __init__(self, filename: str):
        self._filename: str = filename
        self._mapping = load_file_mappings(filename)


    def mapping(self):
        return self._mapping


def load_file_mappings(filename: str) -> Dict[str, str]:
    with open(filename, mode='rt', encoding='utf-8') as file:
        csv.reader(file, delimiter='\t')

        return []