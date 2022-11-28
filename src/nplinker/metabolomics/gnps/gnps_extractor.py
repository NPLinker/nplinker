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
        self._extract_file_mappings()

    def _extract_molecular_families(self):       
        prefix = "networkedges_selfloop"
        utils.extract_file_matching_pattern(
            self.data(),
            prefix,
            ".pairsinfo",
            self.target(),
            "molecular_families.pairsinfo"
        )
        os.rmdir(self.target() / prefix)
    
    def _extract_file_mappings(self):
        prefix = "clusterinfosummarygroup_attributes_withIDs_withcomponentID"
        suffix = ".tsv"
        utils.extract_file_matching_pattern(
            self.data(),
            prefix,
            suffix,
            self.target(),
            "file_mappings.tsv"
        )
        os.rmdir(self.target() / prefix)
    

        