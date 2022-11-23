import os
from pathlib import Path
import zipfile


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
        spectra_out = self.target() / "spectra.mgf"
        
        files: list[Path] = [x.filename for x in self.data().filelist]
        spectra_file = list(filter(lambda x: x.endswith(".mgf"), files))[0]

        self.data().extract(spectra_file, self.target())
        os.rename(self.target() / spectra_file, spectra_out)

        