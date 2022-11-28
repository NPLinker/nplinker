import os
from pathlib import Path
import zipfile

from nplinker import utils
from nplinker.metabolomics.gnps.gnps_format import GNPSFormat, gnps_format_from_archive


class GNPSExtractor:
    def __init__(self, filepath: Path, extract_path: Path):
        self._filepath = filepath
        self._extract_path = extract_path
        self._is_fbmn = gnps_format_from_archive(self.data()) == GNPSFormat.FBMN


    def data(self) -> zipfile.ZipFile:
        """Return the managed archive.

        Returns:
            zipfile.ZipFile: Archive from which data is loaded.
        """
        return zipfile.ZipFile(self._filepath)

    def target(self) -> Path:
        return self._extract_path

    def extract(self):
        self._extract_spectra()
        self._extract_molecular_families()
        self._extract_file_mappings()

    def _extract_spectra(self):
        prefix = "spectra" if self._is_fbmn else ""            
        utils.extract_file_matching_pattern(self.data(), prefix, ".mgf", self.target(), "spectra.mgf")
        if self._is_fbmn:
            os.rmdir(self.target() / prefix)

    def _extract_molecular_families(self):       
        prefix = "networkedges_selfloop"
        suffix = "..selfloop" if self._is_fbmn else ".pairsinfo"
        utils.extract_file_matching_pattern(
            self.data(),
            prefix,
            suffix,
            self.target(),
            "molecular_families.pairsinfo"
        )
        os.rmdir(self.target() / prefix)
    
    def _extract_file_mappings(self):
        prefix = "quantification_table_reformatted" if self._is_fbmn else "clusterinfosummarygroup_attributes_withIDs_withcomponentID"
        suffix = ".csv" if self._is_fbmn else ".tsv"
        utils.extract_file_matching_pattern(
            self.data(),
            prefix,
            suffix,
            self.target(),
            "file_mappings" + suffix
        )
        os.rmdir(self.target() / prefix)
    

        