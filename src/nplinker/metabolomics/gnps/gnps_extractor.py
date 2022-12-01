from os import PathLike
import os
from pathlib import Path
import zipfile

from nplinker import utils
from nplinker.metabolomics.gnps.gnps_format import GNPSFormat, gnps_format_from_archive


class GNPSExtractor:
    def __init__(self, filepath: str | PathLike, extract_path: str | PathLike):
        self._filepath: Path = Path(filepath)
        self._extract_path: Path = Path(extract_path)
        self._is_fbmn = gnps_format_from_archive(self.data()) == GNPSFormat.FBMN


    def data(self) -> zipfile.ZipFile:
        """Return the managed archive.

        Returns:
            zipfile.ZipFile: Archive from which data is loaded.
        """
        return zipfile.ZipFile(self._filepath)

    def target(self) -> str:
        return str(self._extract_path)

    def extract(self):
        self._extract_spectra()
        self._extract_molecular_families()
        self._extract_file_mappings()

    def _extract_spectra(self):
        prefix = "spectra" if self._is_fbmn else ""            
        utils.extract_file_matching_pattern(self.data(), prefix, ".mgf", self._extract_path, "spectra.mgf")
        if self._is_fbmn:
            os.rmdir(self._extract_path / prefix)

    def _extract_molecular_families(self):       
        prefix = "networkedges_selfloop"
        suffix = "..selfloop" if self._is_fbmn else ".pairsinfo"
        utils.extract_file_matching_pattern(
            self.data(),
            prefix,
            suffix,
            self._extract_path,
            "molecular_families.pairsinfo"
        )
        os.rmdir(self._extract_path / prefix)
    
    def _extract_file_mappings(self):
        prefix = "quantification_table_reformatted" if self._is_fbmn else "clusterinfosummarygroup_attributes_withIDs_withcomponentID"
        suffix = ".csv" if self._is_fbmn else ".tsv"
        utils.extract_file_matching_pattern(
            self.data(),
            prefix,
            suffix,
            self._extract_path,
            "file_mappings" + suffix
        )
        os.rmdir(self._extract_path / prefix)
    

        