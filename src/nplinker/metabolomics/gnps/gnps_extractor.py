import os
from pathlib import Path
import zipfile

from nplinker import utils


class GNPSExtractor:
    def __init__(self, filepath: Path, extract_path: Path):
        self._filepath = filepath
        self._extract_path = extract_path
    
    def data(self) -> zipfile.ZipFile:
        """Return the managed archive.

        Returns:
            zipfile.ZipFile: Archive from which data is loaded.
        """
        return zipfile.ZipFile(self._filepath)

    def target(self) -> Path:
        return self._extract_path

    def extract(self):
        utils.extract_file_matching_pattern(self.data(), "", ".mgf", self.target(), "spectra.mgf")
        self._extract_molecular_families()

    def _extract_molecular_families(self):       
        prefix = "networkedges_selfloop"
        utils.extract_file_matching_pattern(self.data(), prefix, ".pairsinfo", self.target(), "molecular_families.pairsinfo")
        os.rmdir(self.target() / prefix)

    

        