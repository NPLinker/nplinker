import fnmatch
import json
import shutil
from glob import glob
from pathlib import Path
import nplinker.globals as globals
from nplinker.config import config
from nplinker.genomics.mibig import download_and_extract_mibig_metadata
from nplinker.metabolomics.gnps import GNPSDownloader
from nplinker.metabolomics.gnps import GNPSExtractor
from nplinker.pairedomics import podp_download_and_extract_antismash_data
from nplinker.pairedomics.runbigscape import run_bigscape
from nplinker.schemas import validate_podp_json
from nplinker.utils import download_url
from nplinker.utils import list_dirs
from nplinker.utils import list_files


PODP_PROJECT_URL = "https://pairedomicsdata.bioinformatics.nl/api/projects/{}"


class DatasetArranger:
    def __init__(self) -> None:
        """Arrange the dataset required by NPLinker.

        This class is used to arrange the datasets required by NPLinker according to the
        configuration. The datasets include MIBiG, GNPS, antiSMASH, and BiG-SCAPE.

        If `config.mode` is "local", the datasets are validated.
        If `config.mode` is "podp", the datasets are downloaded or generated.

        It uses the default downloads directory `globals.DOWNLOADS_DEFAULT_PATH` to store the
        downloaded files. Default data paths for MIBiG, GNPS, antiSMASH, and BiG-SCAPE are defined
        in `nplinker.globals`.
        """
        # Prepare the downloads directory and/or PODP json file which are required for other methods
        globals.DOWNLOADS_DEFAULT_PATH.mkdir(exist_ok=True)
        self.arrange_podp_project_json()

    def arrange(self) -> None:
        """Arrange the datasets according to the configuration.

        The datasets include MIBiG, GNPS, antiSMASH, and BiG-SCAPE.
        """
        self.arrange_mibig()
        self.arrange_gnps()
        self.arrange_antismash()
        self.arrange_bigscape()

    def arrange_podp_project_json(self) -> None:
        """Arrange the PODP project JSON file.

        If `config.mode` is "podp", download the PODP project JSON file if it doesn't exist. Then
        validate the PODP project JSON file if it exists or is downloaded.

        The validation is controlled by the json schema `schemas/podp_adapted_schema.json`.
        """
        if config.mode == "podp":
            podp_file = globals.DOWNLOADS_DEFAULT_PATH / f"paired_datarecord_{config.podp_id}.json"
            if not podp_file.exists():
                download_url(
                    PODP_PROJECT_URL.format(config.podp_id), globals.DOWNLOADS_DEFAULT_PATH
                )

            with open(podp_file, "r") as f:
                json_data = json.load(f)
            validate_podp_json(json_data)

    def arrange_mibig(self) -> None:
        """Arrange the MIBiG metadata.

        Always download and extract the MIBiG metadata if `config.mibig.to_use` is True.
        If the default directory has already existed, it will be removed and re-downloaded to ensure
        the latest version is used. So it's not allowed to manually put MIBiG metadata in the
        default directory.
        """
        if config.mibig.to_use:
            if globals.MIBIG_DEFAULT_PATH.exists():
                # remove existing mibig data
                shutil.rmtree(globals.MIBIG_DEFAULT_PATH)
            download_and_extract_mibig_metadata(
                globals.DOWNLOADS_DEFAULT_PATH,
                globals.MIBIG_DEFAULT_PATH,
                version=config.mibig.version,
            )

    def arrange_gnps(self) -> None:
        """Arrange the GNPS data.

        If `config.mode` is "local", validate the GNPS data directory.
        If `config.mode` is "podp", download the GNPS data if it doesn't exist or remove the
        existing GNPS data and re-download it if it is invalid.

        The validation process includes:
        - Check if the GNPS data directory exists.
        - Check if the required files exist in the GNPS data directory, including:
            - file_mappings.tsv or file_mappings.csv
            - spectra.mgf
            - molecular_families.tsv
            - annotations.tsv
        """
        if config.mode == "local":
            validate_gnps(globals.GNPS_DEFAULT_PATH)

        if config.mode == "podp":
            # set range 3 to ensure download can try 2 times and downloaded data is valid
            for _ in range(3):
                try:
                    validate_gnps(globals.GNPS_DEFAULT_PATH)
                except (FileNotFoundError, ValueError):
                    # Don't need to remove downloaded archive, as it'll be overwritten
                    shutil.rmtree(globals.GNPS_DEFAULT_PATH, ignore_errors=True)
                    self._download_and_extract_gnps()

    def _download_and_extract_gnps(self) -> None:
        """Download and extract the GNPS data.

        Get the GNPS task ID from the PODP project JSON file, then download and extract the GNPS
        data to the default GNPS directory.
        """
        podp_file = globals.DOWNLOADS_DEFAULT_PATH / f"paired_datarecord_{config.podp_id}.json"
        with open(podp_file, "r") as f:
            podp_json_data = json.load(f)
        gnps_task_id = podp_json_data["metabolomics"]["project"].get("molecular_network")

        data_archive = (
            GNPSDownloader(gnps_task_id, globals.DOWNLOADS_DEFAULT_PATH)
            .download()
            .get_download_file()
        )
        GNPSExtractor(data_archive, globals.GNPS_DEFAULT_PATH)

    def arrange_antismash(self) -> None:
        """Arrange the antiSMASH data.

        If `config.mode` is "local", validate the antiSMASH data directory.
        If `config.mode` is "podp", download the antiSMASH data if it doesn't exist or remove the
        existing antiSMASH data and re-download it if it is invalid.

        The validation process includes:
        - Check if the antiSMASH data directory exists.
        - Check if the antiSMASH data directory contains at least one sub-directory, and each
            sub-directory contains at least one BGC file (with the suffix ".region???.gbk" where ???
            is a number).

        AntiSMASH BGC directory must follow the structure below:
        antismash
            ├── genome_id_1 (one AntiSMASH output, e.g. GCF_000514775.1)
            │  ├── GCF_000514775.1.gbk
            │  ├── NZ_AZWO01000004.region001.gbk
            │  └── ...
            ├── genome_id_2
            │  ├── ...
            └── ...
        """
        if config.mode == "local":
            validate_antismash(globals.ANTISMASH_DEFAULT_PATH)

        if config.mode == "podp":
            # set range 3 to ensure download can try 2 times and downloaded data is valid
            for _ in range(3):
                try:
                    validate_antismash(globals.ANTISMASH_DEFAULT_PATH)
                    break
                except FileNotFoundError:
                    shutil.rmtree(globals.ANTISMASH_DEFAULT_PATH, ignore_errors=True)
                    self._download_and_extract_antismash()

    def _download_and_extract_antismash(self) -> None:
        """Download and extract the antiSMASH data.

        Get the antiSMASH data from the PODP project JSON file, then download and extract the
        antiSMASH data to the default antiSMASH directory.
        """
        podp_file = globals.DOWNLOADS_DEFAULT_PATH / f"paired_datarecord_{config.podp_id}.json"
        with open(podp_file, "r") as f:
            podp_json_data = json.load(f)
        podp_download_and_extract_antismash_data(
            podp_json_data["genomes"], globals.DOWNLOADS_DEFAULT_PATH, config.root_dir
        )

    def arrange_bigscape(self) -> None:
        """Arrange the BiG-SCAPE data.

        If `config.mode` is "local", validate the BiG-SCAPE data directory.
        If `config.mode` is "podp", run BiG-SCAPE to generate the clustering file if it doesn't
        exist or remove the existing BiG-SCAPE data and re-run BiG-SCAPE if it is invalid.
        The running output of BiG-SCAPE will be saved to the directory "bigscape_running_output"
        in the default BiG-SCAPE directory, and the clustering file "mix_clustering_c{config.bigscape.cutoff}.tsv"
        will be copied to the default BiG-SCAPE directory.

        The validation process includes:
        - Check if the default BiG-SCAPE data directory exists.
        - Check if the clustering file "mix_clustering_c{config.bigscape.cutoff}.tsv" exists in the
            BiG-SCAPE data directory.
        """
        if config.mode == "local":
            validate_bigscape(globals.BIGSCAPE_DEFAULT_PATH)

        if config.mode == "podp":
            for _ in range(3):
                try:
                    validate_bigscape(globals.BIGSCAPE_DEFAULT_PATH)
                    break
                except FileNotFoundError:
                    shutil.rmtree(globals.BIGSCAPE_DEFAULT_PATH, ignore_errors=True)
                    self._run_bigscape()

    def _run_bigscape(self) -> None:
        """Run BiG-SCAPE to generate the clustering file.

        The output of running BiG-SCAPE will be saved to the directory "bigscape_running_output" in
        the default BiG-SCAPE directory.
        The clustering file "mix_clustering_c{config.bigscape.cutoff}.tsv" will be copied to the
        default BiG-SCAPE directory.
        """
        output = globals.BIGSCAPE_DEFAULT_PATH / "bigscape_running_output"
        output.mkdir(exist_ok=True, parents=True)
        run_bigscape(
            globals.ANTISMASH_DEFAULT_PATH,
            output,
            config.bigscape.parameters,
        )
        for f in glob(str(output / "network_files" / "*" / "mix" / "mix_clustering_c*.tsv")):
            shutil.copy(f, globals.BIGSCAPE_DEFAULT_PATH)


def validate_gnps(gnps_dir: Path) -> None:
    """Validate the GNPS data directory and its contents.

    The GNPS data directory must contain the following files:
    - file_mappings.tsv or file_mappings.csv
    - spectra.mgf
    - molecular_families.tsv
    - annotations.tsv

    Args:
        gnps_dir (Path): Path to the GNPS data directory.

    Raises:
        FileNotFoundError: If the GNPS data directory is not found or any of the required files
            is not found.
        ValueError: If both file_mappings.tsv and file_mapping.csv are found.
    """
    if not gnps_dir.exists():
        raise FileNotFoundError(f"GNPS data directory not found at {gnps_dir}")

    file_mappings_tsv = gnps_dir / "file_mappings.tsv"
    file_mappings_csv = gnps_dir / "file_mappings.csv"
    if file_mappings_tsv.exists() and file_mappings_csv.exists():
        raise ValueError(
            f"Both file_mappings.tsv and file_mapping.csv found in GNPS directory {gnps_dir},"
            f" only one is allowed."
        )
    elif not file_mappings_tsv.exists() and not file_mappings_csv.exists():
        raise FileNotFoundError(
            f"Neither file_mappings.tsv nor file_mapping.csv found in GNPS directory {gnps_dir}"
        )

    required_files = [
        gnps_dir / "spectra.mgf",
        gnps_dir / "molecular_families.tsv",
        gnps_dir / "annotations.tsv",
    ]
    list_not_found = [f.name for f in required_files if not f.exists()]
    if list_not_found:
        raise FileNotFoundError(
            f"Files not found in GNPS directory {gnps_dir}: ', '.join({list_not_found})"
        )


def validate_antismash(antismash_dir: Path) -> None:
    """Validate the antiSMASH data directory and its contents.

    The validation only checks the structure of the antiSMASH data directory and file names.
    It does not check
    - the content of the BGC files
    - the consistency between the antiSMASH data and the PODP project JSON file for the PODP
        mode

    The antiSMASH data directory must exist and contain at least one sub-directory. The name of the
    sub-directories must not contain any space. Each sub-directory must contain at least one BGC
    file (with the suffix ".region???.gbk" where ??? is the region number).

    Args:
        antismash_dir (Path): Path to the antiSMASH data directory.

    Raises:
        FileNotFoundError: If the antiSMASH data directory is not found, or no sub-directories
            are found in the antiSMASH data directory, or no BGC files are found in any
            sub-directory.
        ValueError: If any sub-directory name contains a space.
    """
    if not antismash_dir.exists():
        raise FileNotFoundError(f"antiSMASH data directory not found at {antismash_dir}")

    sub_dirs = list_dirs(antismash_dir)
    if not sub_dirs:
        raise FileNotFoundError(
            "No BGC directories found in antiSMASH data directory {antismash_dir}"
        )

    for sub_dir in sub_dirs:
        dir_name = Path(sub_dir).name
        if " " in dir_name:
            raise ValueError(
                f"antiSMASH sub-directory name {dir_name} contains space, which is not allowed"
            )

        gbk_files = list_files(sub_dir, suffix=".gbk", keep_parent=False)
        bgc_files = fnmatch.filter(gbk_files, "*.region???.gbk")
        if not bgc_files:
            raise FileNotFoundError(f"No BGC files found in antiSMASH sub-directory {sub_dir}")


def validate_bigscape(bigscape_dir: Path) -> None:
    """Validate the BiG-SCAPE data directory and its contents.

    The BiG-SCAPE data directory must exist and contain the clustering file
    "mix_clustering_c{config.bigscape.cutoff}.tsv" where {config.bigscape.cutoff} is the
    bigscape cutoff value set in the config file.

    Args:
        bigscape_dir(Path): Path to the BiG-SCAPE data directory.

    Raises:
        FileNotFoundError: If the BiG-SCAPE data directory or the clustering file is not found.
    """
    if not bigscape_dir.exists():
        raise FileNotFoundError(f"BiG-SCAPE data directory not found at {bigscape_dir}")

    clustering_file = bigscape_dir / f"mix_clustering_c{config.bigscape.cutoff}.tsv"
    if not clustering_file.exists():
        raise FileNotFoundError(f"BiG-SCAPE clustering file not found: {clustering_file}")
