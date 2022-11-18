from pathlib import Path
import httpx


class GNPSDownloader:
    GNPS_DATA_DOWNLOAD_URL = 'https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task={}&view=download_clustered_spectra'

    def __init__(self, task_id: str, outpath: Path):
        """Init the class with a gnps_task_id and a output path on where to save the data.

        Args:
            task_id(str): GNPS task id, identifying the data to be downloaded.
            outpath(Path): Path where to store the downloaded archive.

        Examples:
            >>> GNPSDownloader("c22f44b14a3d450eb836d607cb9521bb", "archive.zip")
            """
        self._task_id = task_id
        self._outpath = outpath

    def download(self):
        """Execute the downloading process. """
        with open(self._outpath, 'wb') as f:
            with httpx.stream('POST', self.url()) as r:
                for data in r.iter_bytes():
                    f.write(data)
   
    
    def task_id(self) -> str:
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
        return GNPSDownloader.GNPS_DATA_DOWNLOAD_URL.format(self._task_id)
    
    