from pathlib import Path
import zipfile


class GNPSLoader:
    def __init__(self, filepath: Path):
        self._filepath = filepath
    
    def data(self) -> zipfile.ZipFile:
        """Return the managed archive.

        Returns:
            zipfile.ZipFile: Archive from which data is loaded.
        """
        return zipfile.ZipFile(self._filepath)