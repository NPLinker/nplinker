from os import PathLike
import os
from pathlib import Path
import zipfile

from nplinker import utils
from nplinker.metabolomics.gnps.gnps_format import GNPSFormat, gnps_format_from_archive


class GNPSExtractor:
    def __init__(self, file: str | PathLike, extract_path: str | PathLike):
        """Class to extract files from a GNPS output archive.
        
        The GNPS output archive contains files for spectra, molecular families and file mappings.

        Args:
            file(str | PathLike): str or PathLike object pointing to the GNPS archive.
            extract_path(str | PathLike): str or PathLike object pointing to where to extract the files to.
        """
        self._file: Path = Path(file)
        self._extract_path: Path = Path(extract_path)
        self._is_fbmn = gnps_format_from_archive(self.data()) == GNPSFormat.FBMN


    def data(self) -> zipfile.ZipFile:
        """Return the managed archive.

        Returns:
            zipfile.ZipFile: Archive from which data is loaded.
        """
        return zipfile.ZipFile(self._file)

    def get_extract_path(self) -> str:
        """Get the path where to extract the files to.

        Returns:
            str: Path where to extract files as string.
        """
        return str(self._extract_path)

    def extract(self):
        """Extract the spectra, molecular family and file mappings file from the handled archive. """
        self._extract_spectra()
        self._extract_molecular_families()
        self._extract_file_mappings()
        self._extract_annotations()

    def _extract_spectra(self):
        """Helper function to extract the spectra file from the archive."""
        prefix = "spectra" if self._is_fbmn else ""            
        utils.extract_file_matching_pattern(self.data(), prefix, ".mgf", self._extract_path, "spectra.mgf")
        if self._is_fbmn:
            os.rmdir(self._extract_path / prefix)

    def _extract_molecular_families(self):
        """Helper function to extract the molecular families file from the archive."""
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
        """Helper function to extract the file mappings file from the archive."""
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

    def _extract_annotations(self):
        """ Helper function to extract the annotations file from the archive. """
        prefix = "DB_result" if self._is_fbmn else "result_specnets_DB"            
        utils.extract_file_matching_pattern(
            self.data(),
            prefix,
            ".tsv",
            self._extract_path,
            "annotations.tsv"
        )
        if self._is_fbmn:
            os.rmdir(self._extract_path / prefix)
    

        