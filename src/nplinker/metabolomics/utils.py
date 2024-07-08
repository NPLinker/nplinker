from __future__ import annotations
import json
import logging
from collections.abc import Mapping
from collections.abc import Sequence
from os import PathLike
from pathlib import Path
from nplinker.schemas import validate_podp_json
from nplinker.strain import StrainCollection
from .gnps.gnps_file_mapping_loader import GNPSFileMappingLoader
from .molecular_family import MolecularFamily
from .spectrum import Spectrum


logger = logging.getLogger(__name__)


def add_annotation_to_spectrum(
    annotations: Mapping[str, dict], spectra: Sequence[Spectrum]
) -> None:
    """Add annotations to the `Spectrum.gnps_annotations` attribute for input spectra.

    It is possible that some spectra don't have annotations.

    !!! note
        The input `spectra` list is changed in place.

    Args:
        annotations: A dictionary of GNPS annotations, where the keys are
            spectrum ids and the values are GNPS annotations.
        spectra: A list of Spectrum objects.
    """
    for spec in spectra:
        if spec.id in annotations:
            spec.gnps_annotations = annotations[spec.id]


def add_strains_to_spectrum(
    strains: StrainCollection, spectra: Sequence[Spectrum]
) -> tuple[list[Spectrum], list[Spectrum]]:
    """Add `Strain` objects to the `Spectrum.strains` attribute for input spectra.

    !!! note
        The input `spectra` list is changed in place.

    Args:
        strains: A collection of strain objects.
        spectra: A list of Spectrum objects.

    Returns:
        A tuple of two lists of Spectrum objects,

            - the first list contains Spectrum objects that are updated with Strain objects;
            - the second list contains Spectrum objects that are not updated with Strain objects
            because no Strain objects are found.
    """
    spectra_with_strains = []
    spectra_without_strains = []
    for spec in spectra:
        try:
            strain_list = strains.lookup(spec.id)
        except ValueError:
            spectra_without_strains.append(spec)
            continue

        for strain in strain_list:
            spec.strains.add(strain)
        spectra_with_strains.append(spec)

    logger.info(
        f"{len(spectra_with_strains)} Spectrum objects updated with Strain objects.\n"
        f"{len(spectra_without_strains)} Spectrum objects not updated with Strain objects."
    )

    return spectra_with_strains, spectra_without_strains


def add_spectrum_to_mf(
    spectra: Sequence[Spectrum], mfs: Sequence[MolecularFamily]
) -> tuple[list[MolecularFamily], list[MolecularFamily], dict[MolecularFamily, set[str]]]:
    """Add Spectrum objects to MolecularFamily objects.

    The attribute `MolecularFamily.spectra_ids` contains the ids of `Spectrum` objects.
    These ids are used to find `Spectrum` objects from the input `spectra` list. The found `Spectrum`
    objects are added to the `MolecularFamily.spectra` attribute.

    It is possible that some spectrum ids are not found in the input `spectra` list, and so their
    `Spectrum` objects are missing in the `MolecularFamily` object.


    !!! note
        The input `mfs` list is changed in place.

    Args:
        spectra: A list of Spectrum objects.
        mfs: A list of MolecularFamily objects.

    Returns:
        A tuple of three elements,

            - the first list contains `MolecularFamily` objects that are updated with `Spectrum` objects
            - the second list contains `MolecularFamily` objects that are not updated with `Spectrum`
            objects (all `Spectrum` objects are missing).
            - the third is a dictionary containing `MolecularFamily` objects as keys and a set of ids
            of missing `Spectrum` objects as values.
    """
    spec_dict = {spec.id: spec for spec in spectra}
    mf_with_spec = []
    mf_without_spec = []
    mf_missing_spec: dict[MolecularFamily, set[str]] = {}
    for mf in mfs:
        for spec_id in mf.spectra_ids:
            try:
                spec = spec_dict[spec_id]
            except KeyError:
                if mf not in mf_missing_spec:
                    mf_missing_spec[mf] = {spec_id}
                else:
                    mf_missing_spec[mf].add(spec_id)
                continue
            mf.add_spectrum(spec)

        if mf.spectra:
            mf_with_spec.append(mf)
        else:
            mf_without_spec.append(mf)

    logger.info(
        f"{len(mf_with_spec)} MolecularFamily objects updated with Spectrum objects.\n"
        f"{len(mf_without_spec)} MolecularFamily objects not updated with Spectrum objects.\n"
        f"{len(mf_missing_spec)} MolecularFamily objects have missing Spectrum objects."
    )
    return mf_with_spec, mf_without_spec, mf_missing_spec


# ------------------------------------------------------------------------------
# Functions to extract mappings for metabolomics side:
# strain_id <-> MS_filename <-> spectrum_id
# ------------------------------------------------------------------------------
def extract_mappings_strain_id_ms_filename(
    podp_project_json_file: str | PathLike,
) -> dict[str, set[str]]:
    """Extract mappings "strain_id <-> MS_filename".

    Args:
        podp_project_json_file: The path to the PODP project JSON file.

    Returns:
        Key is strain id and value is a set of MS filenames.

    Notes:
        The `podp_project_json_file` is the project JSON file downloaded from
        PODP platform. For example, for project MSV000079284, its json file is
        https://pairedomicsdata.bioinformatics.nl/api/projects/4b29ddc3-26d0-40d7-80c5-44fb6631dbf9.4.

    See Also:
        - [podp_generate_strain_mappings][nplinker.strain.utils.podp_generate_strain_mappings]:
            Generate strain mappings JSON file for PODP pipeline.
    """
    mappings_dict: dict[str, set[str]] = {}
    with open(podp_project_json_file, "r") as f:
        json_data = json.load(f)

    validate_podp_json(json_data)

    # Extract mappings strain id <-> metabolomics filename
    for record in json_data["genome_metabolome_links"]:
        strain_id = record["genome_label"]
        # get the actual filename of the mzXML URL
        filename = Path(record["metabolomics_file"]).name
        if strain_id in mappings_dict:
            mappings_dict[strain_id].add(filename)
        else:
            mappings_dict[strain_id] = {filename}
    return mappings_dict


def extract_mappings_ms_filename_spectrum_id(
    gnps_file_mappings_file: str | PathLike,
) -> dict[str, set[str]]:
    """Extract mappings "MS_filename <-> spectrum_id".

    Args:
        gnps_file_mappings_file: The path to the GNPS file mappings file (csv or tsv).

    Returns:
        Key is MS filename and value is a set of spectrum ids.

    Notes:
        The `gnps_file_mappings_file` is downloaded from GNPS website and named as
        [GNPS_FILE_MAPPINGS_TSV][nplinker.defaults.GNPS_FILE_MAPPINGS_TSV] or
        [GNPS_FILE_MAPPINGS_CSV][nplinker.defaults.GNPS_FILE_MAPPINGS_CSV].
        For more details, see [GNPS data][gnps-data].

    See Also:
        - [GNPSFileMappingLoader][nplinker.metabolomics.gnps.gnps_file_mapping_loader.GNPSFileMappingLoader]:
        Load GNPS file mappings file.
        - [podp_generate_strain_mappings][nplinker.strain.utils.podp_generate_strain_mappings]:
            Generate strain mappings JSON file for PODP pipeline.
    """
    loader = GNPSFileMappingLoader(gnps_file_mappings_file)
    return loader.mapping_reversed


def get_mappings_strain_id_spectrum_id(
    mappings_strain_id_ms_filename: Mapping[str, set[str]],
    mappings_ms_filename_spectrum_id: Mapping[str, set[str]],
) -> dict[str, set[str]]:
    """Get mappings "strain_id <-> spectrum_id".

    Args:
        mappings_strain_id_ms_filename: Mappings
            "strain_id <-> MS_filename".
        mappings_ms_filename_spectrum_id: Mappings
            "MS_filename <-> spectrum_id".

    Returns:
        Key is strain id and value is a set of spectrum ids.


    See Also:
        - `extract_mappings_strain_id_ms_filename`: Extract mappings "strain_id <-> MS_filename".
        - `extract_mappings_ms_filename_spectrum_id`: Extract mappings "MS_filename <-> spectrum_id".
        - [podp_generate_strain_mappings][nplinker.strain.utils.podp_generate_strain_mappings]:
            Generate strain mappings JSON file for PODP pipeline.
    """
    mappings_dict = {}
    for strain_id, ms_filenames in mappings_strain_id_ms_filename.items():
        spectrum_ids = set()
        for ms_filename in ms_filenames:
            if (sid := mappings_ms_filename_spectrum_id.get(ms_filename)) is not None:
                spectrum_ids.update(sid)
        if spectrum_ids:
            mappings_dict[strain_id] = spectrum_ids
    return mappings_dict
