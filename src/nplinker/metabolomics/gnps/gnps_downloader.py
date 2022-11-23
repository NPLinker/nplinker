from pathlib import Path
from bs4 import BeautifulSoup, Tag
import httpx
import requests


class GNPSDownloader:
    GNPS_TASK_URL = 'https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task={}'
    GNPS_DATA_DOWNLOAD_URL = 'https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task={}&view=download_clustered_spectra'
    GNPS_DATA_DOWNLOAD_URL_FBMN = 'https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task={}&view=download_cytoscape_data'

    def __init__(self, task_id: str, download_root: Path):
        """Init the class with a gnps_task_id and a output path on where to save the data.

        Args:
            task_id(str): GNPS task id, identifying the data to be downloaded.
            download_root(Path): Path where to store the downloaded archive.

        Examples:
            >>> GNPSDownloader("c22f44b14a3d450eb836d607cb9521bb", "~/downloads")
            """
        self._task_id = task_id
        self._download_root = download_root


    def download(self):
        """Execute the downloading process. """
        with open(self.target(), 'wb') as f:
            with httpx.stream('POST', self.url()) as r:
                for data in r.iter_bytes():
                    f.write(data)
   
    def target(self) -> Path:
        return self._download_root.joinpath(self._task_id + ".zip")
    
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
        
        task_html = requests.get(GNPSDownloader.GNPS_TASK_URL.format(self._task_id))        
        soup = BeautifulSoup(task_html.text)
        tags = soup.find_all('th')
        workflow_tag: Tag = list(filter(lambda x: x.contents == ['Workflow'], tags))[0]
        workflow_format_tag: Tag = workflow_tag.parent.contents[3]
        workflow_format = workflow_format_tag.contents[0].strip()
        
        
        return GNPSDownloader.GNPS_DATA_DOWNLOAD_URL.format(self._task_id)
    
    