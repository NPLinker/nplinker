import csv
from typing import Dict, Literal
from nplinker.metabolomics.IFileMappingLoader import IFileMappingLoader
from nplinker.metabolomics.load_gnps import identify_gnps_format
from nplinker.utils import find_delimiter


class GNPSFileMappingLoader(IFileMappingLoader):

    def __init__(self, filename: str):
        self._filename: str = filename
        self._mapping: Dict[int, str] = {}

    def mapping(self):
        return self._mapping

    def load_mapping_allfiles(self):
        with open(self._filename, mode='rt', encoding='utf-8') as file:

            delimiter = find_delimiter(self._filename)
            reader = csv.reader(file, delimiter=delimiter)
            header: list[str] = next(reader)
            reader = csv.DictReader(file, header, delimiter=delimiter)
 
            for row in reader:
                spectrum_id = int(row["cluster index"])
                occurrences = row["AllFiles"].split("###") # split by '###'
                occurrences.pop() # remove last empty entry

                # separate the scan position from the files
                samples = [x.split(':')[0] for x in occurrences]

                self._mapping[spectrum_id] = samples

    def load_mapping_fbmn(self):
        with open(self._filename, mode='rt', encoding='utf-8') as file:
            reader = csv.reader(file)
            line: list[str] = next(reader)
