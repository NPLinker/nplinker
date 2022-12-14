import csv
from os import PathLike
from pathlib import Path
from typing import Any

from nplinker.metabolomics.abc import AnnotationLoaderBase

GNPS_URL_FORMAT = 'https://metabolomics-usi.ucsd.edu/{}/?usi=mzspec:GNPSLIBRARY:{}'

class GNPSAnnotationLoader(AnnotationLoaderBase):
    def __init__(self, file: str | PathLike):
        """Load annotations from GNPS output file.

        Args:
            file(str | PathLike): The GNPS annotation file.
        """
        self._file = Path(file)
        self._annotations : dict[int, dict] = dict()

        with open(self._file, mode='rt', encoding='UTF-8') as f:
            header = f.readline().split('\t')
            dict_reader = csv.DictReader(f, header, delimiter='\t')
            for row in dict_reader:
                scan_id = int(row.pop('#Scan#'))
                self._annotations[scan_id] = row
                
                # also insert useful URLs
                for t in ['png', 'json', 'svg', 'spectrum']:
                    self._annotations[scan_id][f'{t}_url'] = GNPS_URL_FORMAT.format(t, row['SpectrumID'])

            

    def annotations(self) -> dict[int, dict]:
        """Get annotations.

        Returns:
            dict[int, dict]: Spectra indices are keys and values are the annotations for this spectrum.

        Examples:
            >>> print(loader.annotations()[100])
            """
        return self._annotations