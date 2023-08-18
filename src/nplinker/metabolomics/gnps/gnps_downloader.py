from os import PathLike
from pathlib import Path
import httpx
from typing_extensions import Self
from nplinker.metabolomics.gnps.gnps_format import gnps_format_from_task_id
from nplinker.metabolomics.gnps.gnps_format import GNPSFormat


class GNPSDownloader:
    GNPS_DATA_DOWNLOAD_URL = 'https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task={}&view=download_clustered_spectra'
    GNPS_DATA_DOWNLOAD_URL_FBMN = 'https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task={}&view=download_cytoscape_data'

    def __init__(self, task_id: str, download_root: str | PathLike):
        """Download GNPS zip archive for the given task id.

        Note that only GNPS workflows listed in the GNPSFormat enum are supported.

        Args:
            task_id(str): GNPS task id, identifying the data to be downloaded.
            download_root(Path): Path where to store the downloaded archive.

        Raises:
            ValueError: If the given task id does not correspond to a supported
                GNPS workflow.

        Examples:
            >>> GNPSDownloader("c22f44b14a3d450eb836d607cb9521bb", "~/downloads")
        """
        gnps_format = gnps_format_from_task_id(task_id)
        if gnps_format == GNPSFormat.Unknown:
            raise ValueError(
                f"Unknown workflow type for GNPS task '{task_id}'."
                f"Supported GNPS workflows are described in the GNPSFormat enum, "
                f"including such as 'METABOLOMICS-SNETS', 'METABOLOMICS-SNETS-V2' "
                f"and 'FEATURE-BASED-MOLECULAR-NETWORKING'.")

        self._task_id = task_id
        self._download_root: Path = Path(download_root)
        self.gnps_format = gnps_format
        self._file_name = gnps_format.value + "-" + self._task_id + ".zip"

    def download(self) -> Self:
        """Execute the downloading process. """
        with open(self.get_download_path(), 'wb') as f:
            with httpx.stream('POST', self.get_url()) as r:
                for data in r.iter_bytes():
                    f.write(data)
        return self

    def get_download_path(self) -> str:
        """Get the path where to store the downloaded file.

        Returns:
            str: Download path as string
        """
        return str(self._download_root.joinpath(self._task_id + ".zip"))

    def get_task_id(self) -> str:
        """Get the GNPS task id.

        Returns:
            str: Task id as string.
        """
        return self._task_id

    def get_url(self) -> str:
        """Get the full URL linking to GNPS data to be dowloaded.

        Returns:
            str: URL pointing to the GNPS data to be downloaded.
        """
        gnps_format = gnps_format_from_task_id(self._task_id)

        if gnps_format == GNPSFormat.Unknown:
            raise ValueError(
                f"Unknown workflow type for GNPS task '{self._task_id}'."
                f"Supported GNPS workflows are: 'METABOLOMICS-SNETS', "
                f"'METABOLOMICS-SNETS-V2', 'FEATURE-BASED-MOLECULAR-NETWORKING'"
            )

        if gnps_format == GNPSFormat.FBMN:
            return GNPSDownloader.GNPS_DATA_DOWNLOAD_URL_FBMN.format(
                self._task_id)

        return GNPSDownloader.GNPS_DATA_DOWNLOAD_URL.format(self._task_id)
