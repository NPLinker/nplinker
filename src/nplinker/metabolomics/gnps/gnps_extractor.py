from __future__ import annotations
import os
import zipfile
from os import PathLike
from pathlib import Path
from nplinker import utils
from .gnps_format import GNPSFormat
from .gnps_format import gnps_format_from_archive


class GNPSExtractor:
    """Extract files from a GNPS molecular networking archive (.zip).

    ??? info "Concept"
        [GNPS data][gnps-data]

    Four files are extracted and renamed to the following names:

    - file_mappings(.tsv/.csv)
    - spectra.mgf
    - molecular_families.tsv
    - annotations.tsv

    The files to be extracted are selected based on the GNPS workflow type,
    as described below (in the order of the files above):

    1. METABOLOMICS-SNETS
        - clusterinfosummarygroup_attributes_withIDs_withcomponentID/*.tsv
        - METABOLOMICS-SNETS*.mgf
        - networkedges_selfloop/*.pairsinfo
        - result_specnets_DB/*.tsv
    2. METABOLOMICS-SNETS-V2
        - clusterinfosummarygroup_attributes_withIDs_withcomponentID/*.clustersummary
        - METABOLOMICS-SNETS-V2*.mgf
        - networkedges_selfloop/*.selfloop
        - result_specnets_DB/.tsv
    3. FEATURE-BASED-MOLECULAR-NETWORKING
        - quantification_table*/*.csv
        - spectra/*.mgf
        - networkedges_selfloop/*.selfloop
        - DB_result/*.tsv

    Attributes:
        gnps_format: The GNPS workflow type.
        extract_dir: The path where to extract the files to.
    """

    def __init__(self, file: str | PathLike, extract_dir: str | PathLike):
        """Initialize the GNPSExtractor.

        Args:
            file: The path to the GNPS zip file.
            extract_dir: path to the directory where to extract the files to.

        Raises:
            ValueError: If the given file is an invalid GNPS archive.

        Examples:
            >>> gnps_extractor = GNPSExtractor("path/to/gnps_archive.zip", "path/to/extract_dir")
            >>> gnps_extractor.gnps_format
            <GNPSFormat.SNETS: 'METABOLOMICS-SNETS'>
            >>> gnps_extractor.extract_dir
            'path/to/extract_dir'
        """
        gnps_format = gnps_format_from_archive(file)
        if gnps_format == GNPSFormat.Unknown:
            raise ValueError(
                f"Unknown workflow type for GNPS archive '{file}'."
                f"Supported GNPS workflows are described in the GNPSFormat enum, "
                f"including such as 'METABOLOMICS-SNETS', 'METABOLOMICS-SNETS-V2' "
                f"and 'FEATURE-BASED-MOLECULAR-NETWORKING'."
            )

        self._file = Path(file)
        self._extract_path = Path(extract_dir)
        self._gnps_format = gnps_format
        # the order of filenames matters
        self._target_files = [
            "file_mappings",
            "spectra.mgf",
            "molecular_families.tsv",
            "annotations.tsv",
        ]

        self._extract()

    @property
    def gnps_format(self) -> GNPSFormat:
        """Get the GNPS workflow type.

        Returns:
            GNPS workflow type.
        """
        return self._gnps_format

    @property
    def extract_dir(self) -> str:
        """Get the path where to extract the files to.

        Returns:
            Path where to extract files as string.
        """
        return str(self._extract_path)

    def _extract(self):
        """Extract required files from archive."""
        if self._gnps_format == GNPSFormat.SNETS:
            self._extract_snets()
        elif self._gnps_format == GNPSFormat.SNETSV2:
            self._extract_snetsv2()
        elif self._gnps_format == GNPSFormat.FBMN:
            self._extract_fbmn()

    def _extract_snets(self):
        # the order of members matters
        members = [
            self._select_member(
                "clusterinfosummarygroup_attributes_withIDs_withcomponentID", ".tsv"
            ),
            self._select_member("METABOLOMICS-SNETS", ".mgf"),
            self._select_member("networkedges_selfloop", ".pairsinfo"),
            self._select_member("result_specnets_DB", ".tsv"),
        ]
        utils.extract_archive(self._file, self._extract_path, members)
        # rename the files to the expected names
        # os.renames automatically remove empty directories after renaming
        os.renames(
            self._extract_path / members[0], self._extract_path / (self._target_files[0] + ".tsv")
        )
        for member, fname in zip(members[1:], self._target_files[1:]):
            os.renames(self._extract_path / member, self._extract_path / fname)

    def _extract_snetsv2(self):
        # the order of members matters
        members = [
            self._select_member(
                "clusterinfosummarygroup_attributes_withIDs_withcomponentID", ".clustersummary"
            ),
            self._select_member("METABOLOMICS-SNETS-V2", ".mgf"),
            self._select_member("networkedges_selfloop", ".selfloop"),
            self._select_member("result_specnets_DB", ".tsv"),
        ]
        utils.extract_archive(self._file, self._extract_path, members)
        os.renames(
            self._extract_path / members[0], self._extract_path / (self._target_files[0] + ".tsv")
        )
        for member, fname in zip(members[1:], self._target_files[1:]):
            os.renames(self._extract_path / member, self._extract_path / fname)

    def _extract_fbmn(self):
        # the order of members matters
        members = [
            self._select_member("quantification_table", ".csv"),
            self._select_member("spectra", ".mgf"),
            self._select_member("networkedges_selfloop", ".selfloop"),
            self._select_member("DB_result", ".tsv"),
        ]
        utils.extract_archive(self._file, self._extract_path, members)
        os.renames(
            self._extract_path / members[0], self._extract_path / (self._target_files[0] + ".csv")
        )
        for member, fname in zip(members[1:], self._target_files[1:]):
            os.renames(self._extract_path / member, self._extract_path / fname)

    def _select_member(self, prefix: str, suffix: str) -> str:
        """Helper function to extract files matching a prefix and suffix from the archive."""
        with zipfile.ZipFile(self._file) as zf:
            member_list = [
                member
                for member in zf.namelist()
                if member.startswith(prefix) and member.endswith(suffix)
            ]
            if len(member_list) != 1:
                raise ValueError(
                    f"Expected exactly one file matching pattern '{prefix}*{suffix}'"
                    f"in archive '{self._file}', but found {len(member_list)}."
                )
        return member_list[0]
