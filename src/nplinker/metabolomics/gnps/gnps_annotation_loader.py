import csv
from os import PathLike
from pathlib import Path
from typing import Any

GNPS_URL_FORMAT = 'https://metabolomics-usi.ucsd.edu/{}/?usi=mzspec:GNPSLIBRARY:{}'

class GNPSAnnotationLoader():
    def __init__(self, file: str | PathLike):
        self._file = Path(file)
        self._annotations : dict[int, dict[Any, Any]] = dict()

        with open(self._file, mode='rt', encoding='UTF-8') as f:
            header = f.readline().split('\t')
            dict_reader = csv.DictReader(f, header, delimiter='\t')
            for row in dict_reader:
                scan_id = int(row.pop('#Scan#'))
                self._annotations[scan_id] = row
                
                # also insert useful URLs
                for t in ['png', 'json', 'svg', 'spectrum']:
                    self._annotations[scan_id][f'{t}_url'] = GNPS_URL_FORMAT.format(t, row['SpectrumID'])

            

    def annotations(self) -> dict[int, dict[Any, Any]]:
        return self._annotations