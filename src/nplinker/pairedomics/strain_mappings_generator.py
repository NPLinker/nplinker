import json
import logging
from os import PathLike
from pathlib import Path
from jsonschema import validate
from nplinker.metabolomics.gnps import GNPSFileMappingLoader
from nplinker.schemas import GENOME_BGC_MAPPINGS_SCHEMA
from nplinker.schemas import validate_podp_json
from nplinker.strain_collection import StrainCollection
from nplinker.strains import Strain
from .podp_antismash_downloader import GenomeStatus
from .podp_antismash_downloader import get_best_available_genome_id


logger = logging.getLogger(__name__)

__all__ = [
    "podp_generate_strain_mappings",
    "extract_mappings_strain_id_original_genome_id",
    "extract_mappings_original_genome_id_resolved_genome_id",
    "extract_mappings_resolved_genome_id_bgc_id",
    "get_mappings_strain_id_bgc_id", "extract_mappings_strain_id_ms_filename",
    "extract_mappings_ms_filename_spectrum_id",
    "get_mappings_strain_id_spectrum_id"
]


def podp_generate_strain_mappings(
        podp_project_json_file: str | PathLike,
        genome_status_json_file: str | PathLike,
        genome_bgc_mappings_file: str | PathLike,
        gnps_file_mapping_tsv_file: str | PathLike,
        output_json_file: str | PathLike) -> StrainCollection:
    """Generate strain mappings JSON file for PODP pipeline.

    To get the strain mappings, we need to combine the following mappings:
    - strain_id <-> original_genome_id <-> resolved_genome_id <-> bgc_id
    - strain_id <-> MS_filename <-> spectrum_id

    These mappings are extracted from the following files:
    - "strain_id <-> original_genome_id" is extracted from `podp_project_json_file`.
    - "original_genome_id <-> resolved_genome_id" is extracted from `genome_status_json_file`.
    - "resolved_genome_id <-> bgc_id" is extracted from `genome_bgc_mappings_file`.
    - "strain_id <-> MS_filename" is extracted from `podp_project_json_file`.
    - "MS_filename <-> spectrum_id" is extracted from `gnps_file_mapping_tsv_file`.

    Args:
        podp_project_json_file(str | PathLike): The path to the PODP project
            JSON file.
        genome_status_json_file(str | PathLike): The path to the genome status
            JSON file.
        genome_bgc_mappings_file(str | PathLike): The path to the genome BGC
            mappings JSON file.
        gnps_file_mapping_tsv_file(str | PathLike): The path to the GNPS file
            mapping TSV file.
        output_json_file(str | PathLike): The path to the output JSON file.

    Returns:
        StrainCollection: The strain mappings stored in a StrainCollection object.

    See Also:
        `extract_mappings_strain_id_original_genome_id`: Extract mappings
            "strain_id <-> original_genome_id".
        `extract_mappings_original_genome_id_resolved_genome_id`: Extract mappings
            "original_genome_id <-> resolved_genome_id".
        `extract_mappings_resolved_genome_id_bgc_id`: Extract mappings
            "resolved_genome_id <-> bgc_id".
        `get_mappings_strain_id_bgc_id`: Get mappings "strain_id <-> bgc_id".
        `extract_mappings_strain_id_ms_filename`: Extract mappings
            "strain_id <-> MS_filename".
        `extract_mappings_ms_filename_spectrum_id`: Extract mappings
            "MS_filename <-> spectrum_id".
        `get_mappings_strain_id_spectrum_id`: Get mappings "strain_id <-> spectrum_id".
    """

    # Get mappings strain_id <-> original_geonme_id <-> resolved_genome_id <-> bgc_id
    mappings_strain_id_bgc_id = get_mappings_strain_id_bgc_id(
        extract_mappings_strain_id_original_genome_id(podp_project_json_file),
        extract_mappings_original_genome_id_resolved_genome_id(
            genome_status_json_file),
        extract_mappings_resolved_genome_id_bgc_id(genome_bgc_mappings_file))

    # Get mappings strain_id <-> MS_filename <-> spectrum_id
    mappings_strain_id_spectrum_id = get_mappings_strain_id_spectrum_id(
        extract_mappings_strain_id_ms_filename(podp_project_json_file),
        extract_mappings_ms_filename_spectrum_id(gnps_file_mapping_tsv_file))

    # Get mappings strain_id <-> bgc_id / spectrum_id
    mappings = mappings_strain_id_bgc_id.copy()
    for strain_id, spectrum_ids in mappings_strain_id_spectrum_id.items():
        if strain_id in mappings:
            mappings[strain_id].update(spectrum_ids)
        else:
            mappings[strain_id] = spectrum_ids.copy()

    # Create StrainCollection
    sc = StrainCollection()
    for strain_id, bgc_ids in mappings.items():
        if not sc.has_name(strain_id):
            strain = Strain(strain_id)
            for bgc_id in bgc_ids:
                strain.add_alias(bgc_id)
            sc.add(strain)
        else:
            # strain_list has only one element
            strain_list = sc.lookup(strain_id)
            for bgc_id in bgc_ids:
                strain_list[0].add_alias(bgc_id)

    # Write strain mappings JSON file
    sc.to_json(output_json_file)
    logger.info('Generated strain mappings JSON file: %s', output_json_file)

    return sc


#------------------------------------------------------------------------------
# Functions to extract mappings for genomics side:
# strain_id <-> original_geonme_id <-> resolved_genome_id <-> bgc_id
#------------------------------------------------------------------------------
def extract_mappings_strain_id_original_genome_id(
        podp_project_json_file: str | PathLike) -> dict[str, set[str]]:
    """Extract mappings "strain id <-> original genome id".

    Args:
        podp_project_json_file(str | PathLike): The path to the PODP project
            JSON file.

    Returns:
        dict[str, set[str]]: Key is strain id and value is a set of original genome ids.

    Notes:
        The `podp_project_json_file` is the project JSON file downloaded from
        PODP platform. For example, for project MSV000079284, its json file is
        https://pairedomicsdata.bioinformatics.nl/api/projects/4b29ddc3-26d0-40d7-80c5-44fb6631dbf9.4.
    """
    mappings_dict = {}
    with open(podp_project_json_file, 'r') as f:
        json_data = json.load(f)

    validate_podp_json(json_data)

    for record in json_data['genomes']:
        strain_id = record['genome_label']
        genome_id = get_best_available_genome_id(record['genome_ID'])
        if genome_id is None:
            logger.warning(
                'Failed to extract genome ID from genome with label %s',
                strain_id)
            continue
        if strain_id in mappings_dict:
            mappings_dict[strain_id].add(genome_id)
        else:
            mappings_dict[strain_id] = {genome_id}
    return mappings_dict


def extract_mappings_original_genome_id_resolved_genome_id(
        genome_status_json_file: str | PathLike) -> dict[str, str]:
    """Extract mappings "original_genome_id <-> resolved_genome_id".

    Args:
        genome_status_json_file(str | PathLike): The path to the genome status
            JSON file.

    Returns:
        dict[str, str]: Key is original genome id and value is resolved genome id.

    Notes:
        The `genome_status_json_file` is usually generated by the
        `podp_download_and_extract_antismash_data` function with
        a default file name defined in `nplinker.globals.GENOME_STATUS_FILENAME`.
    """
    gs_mappings_dict = GenomeStatus.read_json(genome_status_json_file)
    return {
        gs.original_id: gs.resolved_refseq_id
        for gs in gs_mappings_dict.values()
    }


def extract_mappings_resolved_genome_id_bgc_id(
        genome_bgc_mappings_file: str | PathLike) -> dict[str, set[str]]:
    """Extract mappings "resolved_genome_id <-> bgc_id".

    Args:
        genome_bgc_mappings_file(str | PathLike): The path to the genome BGC
            mappings JSON file.

    Returns:
        dict[str, set[str]]: Key is resolved genome id and value is a set of BGC ids.

    Notes:
        The `genome_bgc_mappings_file` is usually generated by the
        `generate_mappings_genome_id_bgc_id` function with a default file name
        defined in `nplinker.globals.GENOME_BGC_MAPPINGS_FILENAME`.
    """
    with open(genome_bgc_mappings_file, 'r') as f:
        json_data = json.load(f)

    # validate the JSON data
    validate(json_data, GENOME_BGC_MAPPINGS_SCHEMA)

    return {
        mapping["genome_ID"]: set(mapping["BGC_ID"])
        for mapping in json_data['mappings']
    }


def get_mappings_strain_id_bgc_id(
    mappings_strain_id_original_genome_id: dict[str, set[str]],
    mappings_original_genome_id_resolved_genome_id: dict[str, str],
    mappings_resolved_genome_id_bgc_id: dict[str, set[str]]
) -> dict[str, set[str]]:
    """Get mappings "strain_id <-> bgc_id".

    Args:
        mappings_strain_id_original_genome_id(dict[str, set[str]]): Mappings
            "strain_id <-> original_genome_id".
        mappings_original_genome_id_resolved_genome_id(dict[str, str]): Mappings
            "original_genome_id <-> resolved_genome_id".
        mappings_resolved_genome_id_bgc_id(dict[str, set[str]]): Mappings
            "resolved_genome_id <-> bgc_id".

    Returns:
        dict[str, set[str]]: Key is strain id and value is a set of BGC ids.

    See Also:
        `extract_mappings_strain_id_original_genome_id`: Extract mappings
            "strain_id <-> original_genome_id".
        `extract_mappings_original_genome_id_resolved_genome_id`: Extract mappings
            "original_genome_id <-> resolved_genome_id".
        `extract_mappings_resolved_genome_id_bgc_id`: Extract mappings
            "resolved_genome_id <-> bgc_id".
    """
    mappings_dict = {}
    for strain_id, original_genome_ids in mappings_strain_id_original_genome_id.items(
    ):
        bgc_ids = set()
        for original_genome_id in original_genome_ids:
            resolved_genome_id = mappings_original_genome_id_resolved_genome_id[
                original_genome_id]
            if (bgc_id :=
                    mappings_resolved_genome_id_bgc_id.get(resolved_genome_id)
                ) is not None:
                bgc_ids.update(bgc_id)
        if bgc_ids:
            mappings_dict[strain_id] = bgc_ids
    return mappings_dict


#------------------------------------------------------------------------------
# Functions to extract mappings for metabolomics side:
# strain_id <-> MS_filename <-> spectrum_id
#------------------------------------------------------------------------------
def extract_mappings_strain_id_ms_filename(
        podp_project_json_file: str | PathLike) -> dict[str, set[str]]:
    """Extract mappings "strain_id <-> MS_filename".

    Args:
        podp_project_json_file(str | PathLike): The path to the PODP project
            JSON file.

    Returns:
        dict[str, set[str]]: Key is strain id and value is a set of MS filenames.

    Notes:
        The `podp_project_json_file` is the project JSON file downloaded from
        PODP platform. For example, for project MSV000079284, its json file is
        https://pairedomicsdata.bioinformatics.nl/api/projects/4b29ddc3-26d0-40d7-80c5-44fb6631dbf9.4.
    """
    mappings_dict = {}
    with open(podp_project_json_file, 'r') as f:
        json_data = json.load(f)

    validate_podp_json(json_data)

    # Extract mappings strain id <-> metabolomics filename
    for record in json_data['genome_metabolome_links']:
        strain_id = record['genome_label']
        # get the actual filename of the mzXML URL
        filename = Path(record['metabolomics_file']).name
        if strain_id in mappings_dict:
            mappings_dict[strain_id].add(filename)
        else:
            mappings_dict[strain_id] = {filename}
    return mappings_dict


def extract_mappings_ms_filename_spectrum_id(
        tsv_file: str | PathLike) -> dict[str, set[str]]:
    """Extract mappings "MS_filename <-> spectrum_id".

    Args:
        tsv_file(str | PathLike): The path to the GNPS file mapping TSV file.

    Returns:
        dict[str, set[str]]: Key is MS filename and value is a set of spectrum ids.

    Notes:
        The `tsv_file` is generated by GNPS molecular networking. It's downloaded
        from GNPS website to a file with a default name defined in
        `GNPS_FILE_MAPPINGS_FILENAME`.

    See Also:
        `GNPSFileMappingLoader`: A class to load GNPS file mapping TSV file.
    """
    loader = GNPSFileMappingLoader(tsv_file)
    return loader.mapping_reversed


def get_mappings_strain_id_spectrum_id(
    mappings_strain_id_ms_filename: dict[str, set[str]],
    mappings_ms_filename_spectrum_id: dict[str,
                                           set[str]]) -> dict[str, set[str]]:
    """Get mappings "strain_id <-> spectrum_id".

    Args:
        mappings_strain_id_ms_filename(dict[str, set[str]]): Mappings
            "strain_id <-> MS_filename".
        mappings_ms_filename_spectrum_id(dict[str, set[str]]): Mappings
            "MS_filename <-> spectrum_id".

    Returns:
        dict[str, set[str]]: Key is strain id and value is a set of spectrum ids.


    See Also:
        `extract_mappings_strain_id_ms_filename`: Extract mappings
            "strain_id <-> MS_filename".
        `extract_mappings_ms_filename_spectrum_id`: Extract mappings
            "MS_filename <-> spectrum_id".
    """
    mappings_dict = {}
    for strain_id, ms_filenames in mappings_strain_id_ms_filename.items():
        spectrum_ids = set()
        for ms_filename in ms_filenames:
            if (sid := mappings_ms_filename_spectrum_id.get(ms_filename)
                ) is not None:
                spectrum_ids.update(sid)
        if spectrum_ids:
            mappings_dict[strain_id] = spectrum_ids
    return mappings_dict
