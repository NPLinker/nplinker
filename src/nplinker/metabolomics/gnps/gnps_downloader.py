from os import PathLike
from pathlib import Path
import httpx

from .gnps_format import GNPSFormat
from .gnps_format import gnps_format_from_task_id
from typing_extensions import Self


class GNPSDownloader:
    GNPS_DATA_DOWNLOAD_URL = 'https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task={}&view=download_clustered_spectra'
    GNPS_DATA_DOWNLOAD_URL_FBMN = 'https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task={}&view=download_cytoscape_data'

    def __init__(self, task_id: str, download_root: str | PathLike):
        """Init the class with a gnps_task_id and a output path on where to save the data.

        Args:
            task_id(str): GNPS task id, identifying the data to be downloaded.
            download_root(Path): Path where to store the downloaded archive.

        Examples:
            >>> GNPSDownloader("c22f44b14a3d450eb836d607cb9521bb", "~/downloads")
            """
        self._task_id = task_id
        self._download_root: Path = Path(download_root)


    def download(self) -> Self:
        """Execute the downloading process. """
        with open(self.get_download_path(), 'wb') as f:
            with httpx.stream('POST', self.url()) as r:
                for data in r.iter_bytes():
                    f.write(data)
        return self
   
    def get_download_path(self) -> str:
        """Get the path where to download the files.

        Returns:
            str: Download path as string
        """
        return self._download_root.joinpath(self._task_id + ".zip")
    
    def get_task_id(self) -> str:
        """Get the GNPS task id.

        Returns:
            str: Task id as string.
        """
        return self._task_id
    
    def url(self) -> str:
        """Get the formatted URL pointing to the data to be downloaded for the GNPS task id.

        Returns:
            str: URL pointing to the data.
        """
        
        gnps_format = gnps_format_from_task_id(self._task_id)

        if gnps_format == GNPSFormat.FBMN:
            return GNPSDownloader.GNPS_DATA_DOWNLOAD_URL_FBMN.format(self._task_id)
                
        return GNPSDownloader.GNPS_DATA_DOWNLOAD_URL.format(self._task_id)


    