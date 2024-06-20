from __future__ import annotations
import fnmatch
import json
import shutil
from glob import glob
from os import PathLike
from pathlib import Path
from dynaconf import Dynaconf
from jsonschema import validate
import nplinker.defaults as defaults
from nplinker.genomics.antismash import podp_download_and_extract_antismash_data
from nplinker.genomics.bigscape.runbigscape import run_bigscape
from nplinker.genomics.mibig import download_and_extract_mibig_metadata
from nplinker.genomics.utils import generate_mappings_genome_id_bgc_id
from nplinker.metabolomics.gnps import GNPSDownloader
from nplinker.metabolomics.gnps import GNPSExtractor
from nplinker.schemas import STRAIN_MAPPINGS_SCHEMA
from nplinker.schemas import USER_STRAINS_SCHEMA
from nplinker.schemas import validate_podp_json
from nplinker.strain.utils import podp_generate_strain_mappings
from nplinker.utils import download_url
from nplinker.utils import list_dirs
from nplinker.utils import list_files


PODP_PROJECT_URL = "https://pairedomicsdata.bioinformatics.nl/api/projects/{}"


class DatasetArranger:
    """Arrange the dataset required by NPLinker.

    This class is used to arrange the datasets required by NPLinker according to the
    configuration. The datasets include MIBiG, GNPS, antiSMASH, and BiG-SCAPE.

    If `self.config.mode` is "local", the datasets are validated.
    If `self.config.mode` is "podp", the datasets are downloaded or generated.

    Attributes:
        config: A Dynaconf object that contains the configuration settings. Check `nplinker.config`
            module for more details.
        root_dir: The root directory of the datasets.
        downloads_dir: The directory to store downloaded files.
        mibig_dir: The directory to store MIBiG metadata.
        gnps_dir: The directory to store GNPS data.
        antismash_dir: The directory to store antiSMASH data.
        bigscape_dir: The directory to store BiG-SCAPE data.
        bigscape_running_output_dir: The directory to store the running output of BiG-SCAPE.
    """

    def __init__(self, config: Dynaconf) -> None:
        """Initialize the DatasetArranger.

        Args:
            config: A Dynaconf object that contains the configuration settings. Check `nplinker.config`
                module for more details.
        """
        self.config = config
        self.root_dir = config.root_dir
        self.downloads_dir = self.root_dir / defaults.DOWNLOADS_DIRNAME
        self.downloads_dir.mkdir(exist_ok=True)

        self.mibig_dir = self.root_dir / defaults.MIBIG_DIRNAME
        self.gnps_dir = self.root_dir / defaults.GNPS_DIRNAME
        self.antismash_dir = self.root_dir / defaults.ANTISMASH_DIRNAME
        self.bigscape_dir = self.root_dir / defaults.BIGSCAPE_DIRNAME
        self.bigscape_running_output_dir = (
            self.bigscape_dir / defaults.BIGSCAPE_RUNNING_OUTPUT_DIRNAME
        )

        self.arrange_podp_project_json()

    def arrange(self) -> None:
        """Arrange the datasets according to the configuration.

        The datasets include MIBiG, GNPS, antiSMASH, and BiG-SCAPE.
        """
        # The order of arranging the datasets matters, as some datasets depend on others
        self.arrange_mibig()
        self.arrange_gnps()
        self.arrange_antismash()
        self.arrange_bigscape()
        self.arrange_strain_mappings()
        self.arrange_strains_selected()

    def arrange_podp_project_json(self) -> None:
        """Arrange the PODP project JSON file.

        If `self.config.mode` is "podp", download the PODP project JSON file if it doesn't exist. Then
        validate the PODP project JSON file if it exists or is downloaded.

        The validation is controlled by the json schema `schemas/podp_adapted_schema.json`.
        """
        if self.config.mode == "podp":
            file_name = f"paired_datarecord_{self.config.podp_id}.json"
            podp_file = self.downloads_dir / file_name
            if not podp_file.exists():
                download_url(
                    PODP_PROJECT_URL.format(self.config.podp_id),
                    self.downloads_dir,
                    file_name,
                )

            with open(podp_file, "r") as f:
                json_data = json.load(f)
            validate_podp_json(json_data)

    def arrange_mibig(self) -> None:
        """Arrange the MIBiG metadata.

        Always download and extract the MIBiG metadata if `self.config.mibig.to_use` is True.
        If the default directory has already existed, it will be removed and re-downloaded to ensure
        the latest version is used. So it's not allowed to manually put MIBiG metadata in the
        default directory.
        """
        if self.config.mibig.to_use:
            if self.mibig_dir.exists():
                # remove existing mibig data
                shutil.rmtree(self.mibig_dir)
            download_and_extract_mibig_metadata(
                self.downloads_dir,
                self.mibig_dir,
                version=self.config.mibig.version,
            )

    def arrange_gnps(self) -> None:
        """Arrange the GNPS data.

        If `self.config.mode` is "local", validate the GNPS data directory.
        If `self.config.mode` is "podp", download the GNPS data if it doesn't exist or remove the
        existing GNPS data and re-download it if it is invalid.

        The validation process includes:

        - Check if the GNPS data directory exists.
        - Check if the required files exist in the GNPS data directory, including:
            - file_mappings.tsv or file_mappings.csv
            - spectra.mgf
            - molecular_families.tsv
            - annotations.tsv
        """
        pass_validation = False
        if self.config.mode == "podp":
            # retry downloading at most 3 times if downloaded data has problems
            for _ in range(3):
                try:
                    validate_gnps(self.gnps_dir)
                    pass_validation = True
                    break
                except (FileNotFoundError, ValueError):
                    # Don't need to remove downloaded archive, as it'll be overwritten
                    shutil.rmtree(self.gnps_dir, ignore_errors=True)
                    self._download_and_extract_gnps()

        if not pass_validation:
            validate_gnps(self.gnps_dir)

        # get the path to file_mappings file (csv or tsv)
        self.gnps_file_mappings_file = self._get_gnps_file_mappings_file()

    def _get_gnps_file_mappings_file(self) -> Path:
        """Get the GNPS file mappings file.

        The GNPS file mappings file can be either a TSV file or a CSV file. This method checks if
        the TSV file or the CSV file exists in the default GNPS directory.

        Returns:
            Path to the GNPS file mappings file.
        """
        file_mappings_tsv = self.gnps_dir / defaults.GNPS_FILE_MAPPINGS_TSV
        file_mappings_csv = self.gnps_dir / defaults.GNPS_FILE_MAPPINGS_CSV

        gnps_file_mappings_file = (
            file_mappings_tsv if file_mappings_tsv.exists() else file_mappings_csv
        )

        return gnps_file_mappings_file  # type: ignore

    def _download_and_extract_gnps(self) -> None:
        """Download and extract the GNPS data.

        Get the GNPS task ID from the PODP project JSON file, then download and extract the GNPS
        data to the default GNPS directory.
        """
        podp_file = self.downloads_dir / f"paired_datarecord_{self.config.podp_id}.json"
        with open(podp_file, "r") as f:
            podp_json_data = json.load(f)
        gnps_task_id = podp_json_data["metabolomics"]["project"].get("molecular_network")

        data_archive = (
            GNPSDownloader(gnps_task_id, self.downloads_dir).download().get_download_file()
        )
        GNPSExtractor(data_archive, self.gnps_dir)

    def arrange_antismash(self) -> None:
        """Arrange the antiSMASH data.

        If `self.config.mode` is "local", validate the antiSMASH data directory.
        If `self.config.mode` is "podp", download the antiSMASH data if it doesn't exist or remove the
        existing antiSMASH data and re-download it if it is invalid.

        The validation process includes:
        - Check if the antiSMASH data directory exists.
        - Check if the antiSMASH data directory contains at least one sub-directory, and each
            sub-directory contains at least one BGC file (with the suffix ".region???.gbk" where ???
            is a number).

        AntiSMASH BGC directory must follow the structure below:
        ```
        antismash
            ├── genome_id_1 (one AntiSMASH output, e.g. GCF_000514775.1)
            │  ├── GCF_000514775.1.gbk
            │  ├── NZ_AZWO01000004.region001.gbk
            │  └── ...
            ├── genome_id_2
            │  ├── ...
            └── ...
        ```
        """
        pass_validation = False
        if self.config.mode == "podp":
            for _ in range(3):
                try:
                    validate_antismash(self.antismash_dir)
                    pass_validation = True
                    break
                except FileNotFoundError:
                    shutil.rmtree(self.antismash_dir, ignore_errors=True)
                    self._download_and_extract_antismash()

        if not pass_validation:
            validate_antismash(self.antismash_dir)

    def _download_and_extract_antismash(self) -> None:
        """Download and extract the antiSMASH data.

        Get the antiSMASH data from the PODP project JSON file, then download and extract the
        antiSMASH data to the default antiSMASH directory.
        """
        podp_file = self.downloads_dir / f"paired_datarecord_{self.config.podp_id}.json"
        with open(podp_file, "r") as f:
            podp_json_data = json.load(f)
        podp_download_and_extract_antismash_data(
            podp_json_data["genomes"], self.downloads_dir, self.root_dir
        )

    def arrange_bigscape(self) -> None:
        """Arrange the BiG-SCAPE data.

        If `self.config.mode` is "local", validate the BiG-SCAPE data directory.
        If `self.config.mode` is "podp", run BiG-SCAPE to generate the clustering file if it doesn't
        exist or remove the existing BiG-SCAPE data and re-run BiG-SCAPE if it is invalid.
        The running output of BiG-SCAPE will be saved to the directory "bigscape_running_output"
        in the default BiG-SCAPE directory, and the clustering file
        "mix_clustering_c{self.config.bigscape.cutoff}.tsv" will be copied to the default BiG-SCAPE
        directory.

        The validation process includes:

        - Check if the default BiG-SCAPE data directory exists.
        - Check if the clustering file "mix_clustering_c{self.config.bigscape.cutoff}.tsv" exists in the
                BiG-SCAPE data directory.
        - Check if the 'data_sqlite.db' file exists in the BiG-SCAPE data directory.
        """
        pass_validation = False
        if self.config.mode == "podp":
            for _ in range(3):
                try:
                    validate_bigscape(self.bigscape_dir, self.config.bigscape.cutoff)
                    pass_validation = True
                    break
                except FileNotFoundError:
                    shutil.rmtree(self.bigscape_dir, ignore_errors=True)
                    self._run_bigscape()

        if not pass_validation:
            validate_bigscape(self.bigscape_dir, self.config.bigscape.cutoff)

    def _run_bigscape(self) -> None:
        """Run BiG-SCAPE to generate the clustering file.

        The running output of BiG-SCAPE will be saved to the `BIGSCAPE_RUNNING_OUTPUT_PATH`.

        The clustering file "mix_clustering_c{self.config.bigscape.cutoff}.tsv" will be copied to the
        default BiG-SCAPE directory.
        """
        self.bigscape_running_output_dir.mkdir(exist_ok=True, parents=True)
        run_bigscape(
            self.antismash_dir,
            self.bigscape_running_output_dir,
            self.config.bigscape.parameters,
        )
        for f in glob(
            str(
                self.bigscape_running_output_dir
                / "network_files"
                / "*"
                / "mix"
                / "mix_clustering_c*.tsv"
            )
        ):
            shutil.copy(f, self.bigscape_dir)

    def arrange_strain_mappings(self) -> None:
        """Arrange the strain mappings file.

        If `self.config.mode` is "local", validate the strain mappings file.
        If `self.config.mode` is "podp", always generate the strain mappings file and validate it.

        The validation checks if the strain mappings file exists and if it is a valid JSON file
        according to the schema defined in `schemas/strain_mappings_schema.json`.
        """
        if self.config.mode == "podp":
            self._generate_strain_mappings()

        self._validate_strain_mappings()

    def _validate_strain_mappings(self) -> None:
        """Validate the strain mappings file.

        The validation process includes:

        - Check if the strain mappings file exists.
        - Check if the strain mappings file is a valid JSON file according to the schema defined in
            `schemas/strain_mappings_schema.json`.

        Raises:
            FileNotFoundError: If the strain mappings file is not found.
            ValidationError: If the strain mappings file is not a valid JSON file according to the
                schema.
        """
        strain_mappings_file = self.root_dir / defaults.STRAIN_MAPPINGS_FILENAME

        if not strain_mappings_file.exists():
            raise FileNotFoundError(f"Strain mappings file not found at {strain_mappings_file}")

        with open(strain_mappings_file, "r") as f:
            json_data = json.load(f)
        # Will raise ValidationError if the JSON file is invalid
        validate(instance=json_data, schema=STRAIN_MAPPINGS_SCHEMA)

    def _generate_strain_mappings(self) -> None:
        """Generate the strain mappings file for the PODP mode."""
        podp_json_file = self.downloads_dir / f"paired_datarecord_{self.config.podp_id}.json"
        genome_status_json_file = self.downloads_dir / defaults.GENOME_STATUS_FILENAME
        genome_bgc_mappings_file = self.antismash_dir / defaults.GENOME_BGC_MAPPINGS_FILENAME
        gnps_file_mapping_file = self.gnps_file_mappings_file
        strain_mappings_file = self.root_dir / defaults.STRAIN_MAPPINGS_FILENAME

        # generate the genome_bgc_mappings_file
        generate_mappings_genome_id_bgc_id(self.antismash_dir)
        # generate the strain_mappings_file
        podp_generate_strain_mappings(
            podp_json_file,
            genome_status_json_file,
            genome_bgc_mappings_file,
            gnps_file_mapping_file,
            strain_mappings_file,
        )

    def arrange_strains_selected(self) -> None:
        """Arrange the strains selected file.

        Validate the strains selected file if it exists.
        The validation checks if the strains selected file is a valid JSON file according to the
        schema defined in `schemas/user_strains.json`.
        """
        strains_selected_file = self.root_dir / defaults.STRAINS_SELECTED_FILENAME
        if strains_selected_file.exists():
            with open(strains_selected_file, "r") as f:
                json_data = json.load(f)
            validate(instance=json_data, schema=USER_STRAINS_SCHEMA)


def validate_gnps(gnps_dir: str | PathLike) -> None:
    """Validate the GNPS data directory and its contents.

    The GNPS data directory must contain the following files:

    - file_mappings.tsv or file_mappings.csv
    - spectra.mgf
    - molecular_families.tsv
    - annotations.tsv

    Args:
        gnps_dir: Path to the GNPS data directory.

    Raises:
        FileNotFoundError: If the GNPS data directory is not found or any of the required files
            is not found.
        ValueError: If both file_mappings.tsv and file_mapping.csv are found.
    """
    gnps_dir = Path(gnps_dir)
    if not gnps_dir.exists():
        raise FileNotFoundError(f"GNPS data directory not found at {gnps_dir}")

    file_mappings_tsv = gnps_dir / defaults.GNPS_FILE_MAPPINGS_TSV
    file_mappings_csv = gnps_dir / defaults.GNPS_FILE_MAPPINGS_CSV
    if file_mappings_tsv.exists() and file_mappings_csv.exists():
        raise ValueError(
            f"Both {file_mappings_tsv.name} and {file_mappings_csv.name} found in GNPS directory "
            f"{gnps_dir}, only one is allowed."
        )
    elif not file_mappings_tsv.exists() and not file_mappings_csv.exists():
        raise FileNotFoundError(
            f"Neither {file_mappings_tsv.name} nor {file_mappings_csv.name} found in GNPS directory"
            f" {gnps_dir}"
        )

    required_files = [
        gnps_dir / defaults.GNPS_SPECTRA_FILENAME,
        gnps_dir / defaults.GNPS_MOLECULAR_FAMILY_FILENAME,
        gnps_dir / defaults.GNPS_ANNOTATIONS_FILENAME,
    ]
    list_not_found = [f.name for f in required_files if not f.exists()]
    if list_not_found:
        raise FileNotFoundError(
            f"Files not found in GNPS directory {gnps_dir}: ', '.join({list_not_found})"
        )


def validate_antismash(antismash_dir: str | PathLike) -> None:
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
        antismash_dir: Path to the antiSMASH data directory.

    Raises:
        FileNotFoundError: If the antiSMASH data directory is not found, or no sub-directories
            are found in the antiSMASH data directory, or no BGC files are found in any
            sub-directory.
        ValueError: If any sub-directory name contains a space.
    """
    antismash_dir = Path(antismash_dir)
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


def validate_bigscape(bigscape_dir: str | PathLike, cutoff: str) -> None:
    """Validate the BiG-SCAPE data directory and its contents.

    The BiG-SCAPE data directory must exist and contain the clustering file
    "mix_clustering_c{self.config.bigscape.cutoff}.tsv" where {self.config.bigscape.cutoff} is the
    bigscape cutoff value set in the config file.

    Alternatively, the directory can contain the BiG-SCAPE database file generated by BiG-SCAPE v2.
    At the moment, all the family assignments in the database will be used, so this database should
    contain results from a single run with the desired cutoff.

    Args:
        bigscape_dir: Path to the BiG-SCAPE data directory.
        cutoff: The BiG-SCAPE cutoff value.

    Raises:
        FileNotFoundError: If the BiG-SCAPE data directory or the clustering file is not found.
    """
    bigscape_dir = Path(bigscape_dir)
    if not bigscape_dir.exists():
        raise FileNotFoundError(f"BiG-SCAPE data directory not found at {bigscape_dir}")

    clustering_file = bigscape_dir / f"mix_clustering_c{cutoff}.tsv"
    database_file = bigscape_dir / "data_sqlite.db"
    if not clustering_file.exists() and not database_file.exists():
        raise FileNotFoundError(f"BiG-SCAPE data not found in {clustering_file} or {database_file}")
