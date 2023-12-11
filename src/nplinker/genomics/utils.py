from __future__ import annotations
import json
from os import PathLike
from pathlib import Path
from jsonschema import validate
from nplinker.globals import GENOME_BGC_MAPPINGS_FILENAME
from nplinker.logconfig import LogConfig
from nplinker.schemas import GENOME_BGC_MAPPINGS_SCHEMA
from nplinker.strain_collection import StrainCollection
from nplinker.utils import list_dirs
from nplinker.utils import list_files
from .bgc import BGC
from .gcf import GCF


logger = LogConfig.getLogger(__name__)


def generate_mappings_genome_id_bgc_id(
    bgc_dir: str | PathLike, output_file: str | PathLike | None = None
) -> None:
    """Generate a file that maps genome id to BGC id.

    Note that the `output_file` will be overwritten if it already exists.

    Args:
        bgc_dir(str | PathLike): The directory has one-layer of subfolders and
            each subfolder contains BGC files in `.gbk` format.
            It assumes that
            - the subfolder name is the genome id (e.g. refseq),
            - the BGC file name is the BGC id.
        output_file(str | PathLike | None): The path to the output file. Note
            that the file will be overwritten if it already exists.
            Defaults to None, in which case the output file will be placed in
            the directory `bgc_dir` with a file name defined in global variable
            `GENOME_BGC_MAPPINGS_FILENAME`.
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
    json_data = [{"genome_ID": k, "BGC_ID": v} for k, v in genome_bgc_mappings.items()]
    json_data = {"mappings": json_data, "version": "1.0"}

    # validate json data
    validate(instance=json_data, schema=GENOME_BGC_MAPPINGS_SCHEMA)

    if output_file is None:
        output_file = bgc_dir / GENOME_BGC_MAPPINGS_FILENAME
    with open(output_file, "w") as f:
        json.dump(json_data, f)
    logger.info("Generated genome-BGC mappings file: %s", output_file)


def add_strain_to_bgc(strains: StrainCollection, bgcs: list[BGC]) -> tuple[list[BGC], list[BGC]]:
    """Assign a Strain object to `BGC.strain` for input BGCs.

    BGC id is used to find the corresponding Strain object. It's possible that
    no Strain object is found for a BGC id.

    Note that the input list `bgcs` will be changed in place.

    Args:
        strains(StrainCollection): A collection of all strain objects.
        bgcs(list[BGC]): A list of BGC objects.

    Returns:
        tuple(list[BGC], list[BGC]): A tuple of two lists of BGC objects. The
            first list contains BGC objects that are updated with Strain object;
            the second list contains BGC objects that are not updated with
            Strain object because no Strain object is found.

    Raises:
        ValueError: Multiple strain objects found for a BGC id.
    """
    bgc_with_strain = []
    bgc_without_strain = []
    for bgc in bgcs:
        try:
            strain_list = strains.lookup(bgc.bgc_id)
        except ValueError:
            bgc_without_strain.append(bgc)
            continue
        if len(strain_list) > 1:
            raise ValueError(
                f"Multiple strain objects found for BGC id '{bgc.bgc_id}'."
                f"BGC object accept only one strain."
            )
        bgc.strain = strain_list[0]
        bgc_with_strain.append(bgc)

    logger.info(
        f"{len(bgc_with_strain)} BGC objects updated with Strain object.\n"
        f"{len(bgc_without_strain)} BGC objects not updated with Strain object."
    )
    return bgc_with_strain, bgc_without_strain


def add_bgc_to_gcf(bgcs: list[BGC], gcfs: list[GCF]) -> None:
    """To add BGC objects to GCF object based on GCF's BGC ids.

    This method changes the lists `bgcs` and `gcfs` in place.

    Args:
        bgcs(list[BGC]): A list of BGC objects.
        gcfs(list[GCF]): A list of GCF objects.

    Raises:
        KeyError: BGC id not found in the list of BGC objects.
    """
    bgc_dict = {bgc.bgc_id: bgc for bgc in bgcs}
    for gcf in gcfs:
        for bgc_id in gcf.bgc_ids:
            try:
                bgc = bgc_dict[bgc_id]
            except KeyError as e:
                raise KeyError(
                    f"BGC id '{bgc_id}' from GCF object '{gcf.gcf_id}' "
                    "not found in the list of BGC objects."
                ) from e
            gcf.add_bgc(bgc)


def get_bgcs_from_gcfs(gcfs: list[GCF]) -> list[BGC]:
    """Get all BGC objects from given GCF objects."""
    s = set()
    for gcf in gcfs:
        s |= gcf.bgcs
    return list(s)


def get_strains_from_bgcs(bgcs: list[BGC]) -> StrainCollection:
    """Get all strain objects from given BGC objects."""
    sc = StrainCollection()
    for bgc in bgcs:
        if bgc.strain is not None:
            sc.add(bgc.strain)
        else:
            logger.warning("Strain is None for BGC %s", bgc.bgc_id)
    return sc
