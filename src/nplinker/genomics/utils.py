from __future__ import annotations
import json
import logging
from collections.abc import Mapping
from collections.abc import Sequence
from os import PathLike
from pathlib import Path
from jsonschema import validate
from nplinker.defaults import GENOME_BGC_MAPPINGS_FILENAME
from nplinker.schemas import GENOME_BGC_MAPPINGS_SCHEMA
from nplinker.schemas import validate_podp_json
from nplinker.strain import StrainCollection
from nplinker.utils import list_dirs
from nplinker.utils import list_files
from ..genomics.antismash.podp_antismash_downloader import GenomeStatus
from ..genomics.antismash.podp_antismash_downloader import get_best_available_genome_id
from .bgc import BGC
from .gcf import GCF


logger = logging.getLogger(__name__)


def generate_mappings_genome_id_bgc_id(
    bgc_dir: str | PathLike, output_file: str | PathLike | None = None
) -> None:
    """Generate a file that maps genome id to BGC id.

    The input `bgc_dir` must follow the structure of the `antismash` directory defined in
    [Working Directory Structure][working-directory-structure], e.g.:
    ```shell
    bgc_dir
        ├── genome_id_1
        │  ├── bgc_id_1.gbk
        │  └── ...
        ├── genome_id_2
        │  ├── bgc_id_2.gbk
        │  └── ...
        └── ...
    ```

    Args:
        bgc_dir: The directory has one-layer of subfolders and each subfolder contains BGC files
            in `.gbk` format.

            It assumes that

            - the subfolder name is the genome id (e.g. refseq),
            - the BGC file name is the BGC id.
        output_file: The path to the output file.
            The file will be overwritten if it already exists.

            Defaults to None, in which case the output file will be placed in
            the directory `bgc_dir` with the file name
            [GENOME_BGC_MAPPINGS_FILENAME][nplinker.defaults.GENOME_BGC_MAPPINGS_FILENAME].
    """
    bgc_dir = Path(bgc_dir)
    genome_bgc_mappings = {}

    for subdir in list_dirs(bgc_dir):
        genome_id = Path(subdir).name
        bgc_files = list_files(subdir, suffix=(".gbk"), keep_parent=False)
        bgc_ids = [bgc_id for f in bgc_files if (bgc_id := Path(f).stem) != genome_id]
        if bgc_ids:
            genome_bgc_mappings[genome_id] = bgc_ids
        else:
            logger.warning("No BGC files found in %s", subdir)

    # sort mappings by genome_id and construct json data
    genome_bgc_mappings = dict(sorted(genome_bgc_mappings.items()))
    json_data_mappings = [{"genome_ID": k, "BGC_ID": v} for k, v in genome_bgc_mappings.items()]
    json_data = {"mappings": json_data_mappings, "version": "1.0"}

    # validate json data
    validate(instance=json_data, schema=GENOME_BGC_MAPPINGS_SCHEMA)

    if output_file is None:
        output_file = bgc_dir / GENOME_BGC_MAPPINGS_FILENAME
    with open(output_file, "w") as f:
        json.dump(json_data, f)
    logger.info("Generated genome-BGC mappings file: %s", output_file)


def add_strain_to_bgc(
    strains: StrainCollection, bgcs: Sequence[BGC]
) -> tuple[list[BGC], list[BGC]]:
    """Assign a Strain object to `BGC.strain` for input BGCs.

    BGC id is used to find the corresponding Strain object. It's possible that
    no Strain object is found for a BGC id.

    !!! Note
        The input `bgcs` will be changed in place.

    Args:
        strains: A collection of all strain objects.
        bgcs: A list of BGC objects.

    Returns:
        A tuple of two lists of BGC objects,

            - the first list contains BGC objects that are updated with Strain object;
            - the second list contains BGC objects that are not updated with
                Strain object because no Strain object is found.

    Raises:
        ValueError: Multiple strain objects found for a BGC id.
    """
    bgc_with_strain = []
    bgc_without_strain = []
    for bgc in bgcs:
        try:
            strain_list = strains.lookup(bgc.id)
        except ValueError:
            bgc_without_strain.append(bgc)
            continue
        if len(strain_list) > 1:
            raise ValueError(
                f"Multiple strain objects found for BGC id '{bgc.id}'."
                f"BGC object accept only one strain."
            )
        bgc.strain = strain_list[0]
        bgc_with_strain.append(bgc)

    logger.info(
        f"{len(bgc_with_strain)} BGC objects updated with Strain object.\n"
        f"{len(bgc_without_strain)} BGC objects not updated with Strain object."
    )
    return bgc_with_strain, bgc_without_strain


def add_bgc_to_gcf(
    bgcs: Sequence[BGC], gcfs: Sequence[GCF]
) -> tuple[list[GCF], list[GCF], dict[GCF, set[str]]]:
    """Add BGC objects to GCF object based on GCF's BGC ids.

    The attribute of `GCF.bgc_ids` contains the ids of BGC objects. These ids
    are used to find BGC objects from the input `bgcs` list. The found BGC
    objects are added to the `bgcs` attribute of GCF object. It is possible that
    some BGC ids are not found in the input `bgcs` list, and so their BGC
    objects are missing in the GCF object.

    !!! note
        This method changes the lists `bgcs` and `gcfs` in place.

    Args:
        bgcs: A list of BGC objects.
        gcfs: A list of GCF objects.

    Returns:
        A tuple of two lists and a dictionary,

            - The first list contains GCF objects that are updated with BGC objects;
            - The second list contains GCF objects that are not updated with BGC objects
                because no BGC objects are found;
            - The dictionary contains GCF objects as keys and a set of ids of missing
                BGC objects as values.
    """
    bgc_dict = {bgc.id: bgc for bgc in bgcs}
    gcf_with_bgc = []
    gcf_without_bgc = []
    gcf_missing_bgc: dict[GCF, set[str]] = {}
    for gcf in gcfs:
        for bgc_id in gcf.bgc_ids:
            try:
                bgc = bgc_dict[bgc_id]
            except KeyError:
                if gcf not in gcf_missing_bgc:
                    gcf_missing_bgc[gcf] = {bgc_id}
                else:
                    gcf_missing_bgc[gcf].add(bgc_id)
                continue
            gcf.add_bgc(bgc)

        if gcf.bgcs:
            gcf_with_bgc.append(gcf)
        else:
            gcf_without_bgc.append(gcf)

    logger.info(
        f"{len(gcf_with_bgc)} GCF objects updated with BGC objects.\n"
        f"{len(gcf_without_bgc)} GCF objects not updated with BGC objects.\n"
        f"{len(gcf_missing_bgc)} GCF objects have missing BGC objects."
    )
    return gcf_with_bgc, gcf_without_bgc, gcf_missing_bgc


def get_mibig_from_gcf(gcfs: Sequence[GCF]) -> tuple[list[BGC], StrainCollection]:
    """Get MIBiG BGCs and strains from GCF objects.

    Args:
        gcfs: A list of GCF objects.

    Returns:
        A tuple of two objects,

            - the first is a list of MIBiG BGC objects used in the GCFs;
            - the second is a StrainCollection object that contains all Strain objects used in the
            GCFs.
    """
    mibig_bgcs_in_use = []
    mibig_strains_in_use = StrainCollection()
    for gcf in gcfs:
        for bgc in gcf.bgcs:
            if bgc.is_mibig():
                mibig_bgcs_in_use.append(bgc)
                if bgc.strain is not None:
                    mibig_strains_in_use.add(bgc.strain)
    return mibig_bgcs_in_use, mibig_strains_in_use


# ------------------------------------------------------------------------------
# Functions to extract mappings for genomics side:
# strain_id <-> original_geonme_id <-> resolved_genome_id <-> bgc_id
# ------------------------------------------------------------------------------
def extract_mappings_strain_id_original_genome_id(
    podp_project_json_file: str | PathLike,
) -> dict[str, set[str]]:
    """Extract mappings "strain_id <-> original_genome_id".

    !!! tip
        The `podp_project_json_file` is the JSON file downloaded from PODP platform.

        For example, for PODP project MSV000079284, its JSON file is
        https://pairedomicsdata.bioinformatics.nl/api/projects/4b29ddc3-26d0-40d7-80c5-44fb6631dbf9.4.

    Args:
        podp_project_json_file: The path to the PODP project
            JSON file.

    Returns:
        Key is strain id and value is a set of original genome ids.

    See Also:
        - [podp_generate_strain_mappings][nplinker.strain.utils.podp_generate_strain_mappings]:
            Generate strain mappings JSON file for PODP pipeline.
    """
    mappings_dict: dict[str, set[str]] = {}
    with open(podp_project_json_file, "r") as f:
        json_data = json.load(f)

    validate_podp_json(json_data)

    for record in json_data["genomes"]:
        strain_id = record["genome_label"]
        genome_id = get_best_available_genome_id(record["genome_ID"])
        if genome_id is None:
            logger.warning("Failed to extract genome ID from genome with label %s", strain_id)
            continue
        if strain_id in mappings_dict:
            mappings_dict[strain_id].add(genome_id)
        else:
            mappings_dict[strain_id] = {genome_id}
    return mappings_dict


def extract_mappings_original_genome_id_resolved_genome_id(
    genome_status_json_file: str | PathLike,
) -> dict[str, str]:
    """Extract mappings "original_genome_id <-> resolved_genome_id".

    !!! tip
        The `genome_status_json_file` is generated by the [podp_download_and_extract_antismash_data]
        [nplinker.genomics.antismash.podp_antismash_downloader.podp_download_and_extract_antismash_data]
        function with a default file name [GENOME_STATUS_FILENAME][nplinker.defaults.GENOME_STATUS_FILENAME].

    Args:
        genome_status_json_file: The path to the genome status JSON file.


    Returns:
        Key is original genome id and value is resolved genome id.

    See Also:
        - [podp_generate_strain_mappings][nplinker.strain.utils.podp_generate_strain_mappings]:
            Generate strain mappings JSON file for PODP pipeline.
    """
    gs_mappings_dict = GenomeStatus.read_json(genome_status_json_file)
    return {gs.original_id: gs.resolved_refseq_id for gs in gs_mappings_dict.values()}


def extract_mappings_resolved_genome_id_bgc_id(
    genome_bgc_mappings_file: str | PathLike,
) -> dict[str, set[str]]:
    """Extract mappings "resolved_genome_id <-> bgc_id".

    !!! tip
        The `genome_bgc_mappings_file` is usually generated by the
        [generate_mappings_genome_id_bgc_id][nplinker.genomics.utils.generate_mappings_genome_id_bgc_id]
        function with a default file name [GENOME_BGC_MAPPINGS_FILENAME][nplinker.defaults.GENOME_BGC_MAPPINGS_FILENAME].

    Args:
        genome_bgc_mappings_file: The path to the genome BGC
            mappings JSON file.

    Returns:
        Key is resolved genome id and value is a set of BGC ids.

    See Also:
        - [podp_generate_strain_mappings][nplinker.strain.utils.podp_generate_strain_mappings]:
            Generate strain mappings JSON file for PODP pipeline.
    """
    with open(genome_bgc_mappings_file, "r") as f:
        json_data = json.load(f)

    # validate the JSON data
    validate(json_data, GENOME_BGC_MAPPINGS_SCHEMA)

    return {mapping["genome_ID"]: set(mapping["BGC_ID"]) for mapping in json_data["mappings"]}


def get_mappings_strain_id_bgc_id(
    mappings_strain_id_original_genome_id: Mapping[str, set[str]],
    mappings_original_genome_id_resolved_genome_id: Mapping[str, str],
    mappings_resolved_genome_id_bgc_id: Mapping[str, set[str]],
) -> dict[str, set[str]]:
    """Get mappings "strain_id <-> bgc_id".

    Args:
        mappings_strain_id_original_genome_id: Mappings "strain_id <-> original_genome_id".
        mappings_original_genome_id_resolved_genome_id: Mappings "original_genome_id <-> resolved_genome_id".
        mappings_resolved_genome_id_bgc_id: Mappings "resolved_genome_id <-> bgc_id".

    Returns:
        Key is strain id and value is a set of BGC ids.

    See Also:
        - `extract_mappings_strain_id_original_genome_id`: Extract mappings
            "strain_id <-> original_genome_id".
        - `extract_mappings_original_genome_id_resolved_genome_id`: Extract mappings
            "original_genome_id <-> resolved_genome_id".
        - `extract_mappings_resolved_genome_id_bgc_id`: Extract mappings
            "resolved_genome_id <-> bgc_id".
        - [podp_generate_strain_mappings][nplinker.strain.utils.podp_generate_strain_mappings]:
            Generate strain mappings JSON file for PODP pipeline.
    """
    mappings_dict = {}
    for strain_id, original_genome_ids in mappings_strain_id_original_genome_id.items():
        bgc_ids = set()
        for original_genome_id in original_genome_ids:
            resolved_genome_id = mappings_original_genome_id_resolved_genome_id[original_genome_id]
            if (bgc_id := mappings_resolved_genome_id_bgc_id.get(resolved_genome_id)) is not None:
                bgc_ids.update(bgc_id)
        if bgc_ids:
            mappings_dict[strain_id] = bgc_ids
    return mappings_dict
