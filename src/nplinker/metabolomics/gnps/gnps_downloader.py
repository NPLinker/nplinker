from pathlib import Path
import httpx


class GNPSDownloader:
    GNPS_DATA_DOWNLOAD_URL = 'https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task={}&view=download_clustered_spectra'

    def __init__(self, task_id: str, outpath: Path):
        self._task_id = task_id
        self._outpath = outpath

    def download(self):
        with open(self._outpath, 'wb') as f:
            with httpx.stream('POST', self.url()) as r:
                for data in r.iter_bytes():
                    f.write(data)
   
    
    def task_id(self):
        return self._task_id
    
    def url(self):
        return GNPSDownloader.GNPS_DATA_DOWNLOAD_URL.format(self._task_id)
    
    