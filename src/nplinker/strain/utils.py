from __future__ import annotations
import json
import logging
from os import PathLike
from jsonschema import validate
from nplinker.schemas import USER_STRAINS_SCHEMA
from ..genomics.utils import extract_mappings_original_genome_id_resolved_genome_id
from ..genomics.utils import extract_mappings_resolved_genome_id_bgc_id
from ..genomics.utils import extract_mappings_strain_id_original_genome_id
from ..genomics.utils import get_mappings_strain_id_bgc_id
from ..metabolomics.utils import extract_mappings_ms_filename_spectrum_id
from ..metabolomics.utils import extract_mappings_strain_id_ms_filename
from ..metabolomics.utils import get_mappings_strain_id_spectrum_id
from .strain import Strain
from .strain_collection import StrainCollection


logger = logging.getLogger(__name__)


def load_user_strains(json_file: str | PathLike) -> set[Strain]:
    """Load user specified strains from a JSON file.

    The JSON file will be validated against the schema
    [USER_STRAINS_SCHEMA][nplinker.schemas.USER_STRAINS_SCHEMA]

    The content of the JSON file could be, for example:
    ```
    {"strain_ids": ["strain1", "strain2"]}
    ```

    Args:
        json_file: Path to the JSON file containing user specified strains.

    Returns:
        A set of user specified strains.
    """
    with open(json_file, "r") as f:
        json_data = json.load(f)

    # validate json data
    validate(instance=json_data, schema=USER_STRAINS_SCHEMA)

    strains = set()
    for strain_id in json_data["strain_ids"]:
        strains.add(Strain(strain_id))

    return strains


def podp_generate_strain_mappings(
    podp_project_json_file: str | PathLike,
    genome_status_json_file: str | PathLike,
    genome_bgc_mappings_file: str | PathLike,
    gnps_file_mappings_file: str | PathLike,
    output_json_file: str | PathLike,
) -> StrainCollection:
    """Generate strain mappings JSON file for PODP pipeline.

    To get the strain mappings, we need to combine the following mappings:

    - strain_id <-> original_genome_id <-> resolved_genome_id <-> bgc_id
    - strain_id <-> MS_filename <-> spectrum_id

    These mappings are extracted from the following files:

    - "strain_id <-> original_genome_id" is extracted from `podp_project_json_file`.
    - "original_genome_id <-> resolved_genome_id" is extracted from `genome_status_json_file`.
    - "resolved_genome_id <-> bgc_id" is extracted from `genome_bgc_mappings_file`.
    - "strain_id <-> MS_filename" is extracted from `podp_project_json_file`.
    - "MS_filename <-> spectrum_id" is extracted from `gnps_file_mappings_file`.

    Args:
        podp_project_json_file: The path to the PODP project
            JSON file.
        genome_status_json_file: The path to the genome status
            JSON file.
        genome_bgc_mappings_file: The path to the genome BGC
            mappings JSON file.
        gnps_file_mappings_file: The path to the GNPS file
            mappings file (csv or tsv).
        output_json_file: The path to the output JSON file.

    Returns:
        The strain mappings stored in a StrainCollection object.

    See Also:
        - `extract_mappings_strain_id_original_genome_id`: Extract mappings
            "strain_id <-> original_genome_id".
        - `extract_mappings_original_genome_id_resolved_genome_id`: Extract mappings
            "original_genome_id <-> resolved_genome_id".
        - `extract_mappings_resolved_genome_id_bgc_id`: Extract mappings
            "resolved_genome_id <-> bgc_id".
        - `get_mappings_strain_id_bgc_id`: Get mappings "strain_id <-> bgc_id".
        - `extract_mappings_strain_id_ms_filename`: Extract mappings
            "strain_id <-> MS_filename".
        - `extract_mappings_ms_filename_spectrum_id`: Extract mappings
            "MS_filename <-> spectrum_id".
        - `get_mappings_strain_id_spectrum_id`: Get mappings "strain_id <-> spectrum_id".
    """
    # Get mappings strain_id <-> original_genome_id <-> resolved_genome_id <-> bgc_id
    mappings_strain_id_bgc_id = get_mappings_strain_id_bgc_id(
        extract_mappings_strain_id_original_genome_id(podp_project_json_file),
        extract_mappings_original_genome_id_resolved_genome_id(genome_status_json_file),
        extract_mappings_resolved_genome_id_bgc_id(genome_bgc_mappings_file),
    )

    # Get mappings strain_id <-> MS_filename <-> spectrum_id
    mappings_strain_id_spectrum_id = get_mappings_strain_id_spectrum_id(
        extract_mappings_strain_id_ms_filename(podp_project_json_file),
        extract_mappings_ms_filename_spectrum_id(gnps_file_mappings_file),
    )

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
    logger.info("Generated strain mappings JSON file: %s", output_json_file)

    return sc
