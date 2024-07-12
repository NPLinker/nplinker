from __future__ import annotations
from os import PathLike
from pathlib import Path
from typing_extensions import Self
from nplinker.utils import download_url
from .gnps_format import GNPSFormat
from .gnps_format import gnps_format_from_task_id


class GNPSDownloader:
    """Download GNPS zip archive for the given task id.

    ??? info "Concept"
        [GNPS data][gnps-data]

    Note that only GNPS workflows listed in the GNPSFormat enum are supported.

    Attributes:
        GNPS_DATA_DOWNLOAD_URL: URL template for downloading GNPS data.
        GNPS_DATA_DOWNLOAD_URL_FBMN: URL template for downloading GNPS data for FBMN.
        gnps_format: GNPS workflow type.
    """

    GNPS_DATA_DOWNLOAD_URL: str = (
        "https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task={}&view=download_clustered_spectra"
    )
    GNPS_DATA_DOWNLOAD_URL_FBMN: str = (
        "https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task={}&view=download_cytoscape_data"
    )

    def __init__(self, task_id: str, download_root: str | PathLike):
        """Initialize the GNPSDownloader.

        Args:
            task_id: GNPS task id, identifying the data to be downloaded.
            download_root: Path where to store the downloaded archive.

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
                f"and 'FEATURE-BASED-MOLECULAR-NETWORKING'."
            )

        self._task_id = task_id
        self._download_root: Path = Path(download_root)
        self._gnps_format = gnps_format
        self._file_name = gnps_format.value + "-" + self._task_id + ".zip"

    @property
    def gnps_format(self) -> GNPSFormat:
        """Get the GNPS workflow type.

        Returns:
            GNPS workflow type.
        """
        return self._gnps_format

    def download(self) -> Self:
        """Download GNPS data.

        Note: GNPS data is downloaded using the POST method (empty payload is OK).
        """
        download_url(
            self.get_url(), self._download_root, filename=self._file_name, http_method="POST"
        )
        return self

    def get_download_file(self) -> str:
        """Get the path to the downloaded file.

        Returns:
            Download path as string
        """
        return str(Path(self._download_root) / self._file_name)

    def get_task_id(self) -> str:
        """Get the GNPS task id.

        Returns:
            Task id as string.
        """
        return self._task_id

    def get_url(self) -> str:
        """Get the download URL.

        Returns:
            URL pointing to the GNPS data to be downloaded.
        """
        if self.gnps_format == GNPSFormat.FBMN:
            return GNPSDownloader.GNPS_DATA_DOWNLOAD_URL_FBMN.format(self._task_id)
        return GNPSDownloader.GNPS_DATA_DOWNLOAD_URL.format(self._task_id)
