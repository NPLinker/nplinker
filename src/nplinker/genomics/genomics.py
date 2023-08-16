from __future__ import annotations
import csv
import json
from os import PathLike
from pathlib import Path
from deprecated import deprecated
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
        bgc_dir: str | PathLike,
        output_file: str | PathLike | None = None) -> None:
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
        bgc_ids = [
            bgc_id for f in bgc_files if (bgc_id := Path(f).stem) != genome_id
        ]
        if bgc_ids:
            genome_bgc_mappings[genome_id] = bgc_ids
        else:
            logger.warning("No BGC files found in %s", subdir)

    # sort mappings by genome_id and construct json data
    genome_bgc_mappings = dict(sorted(genome_bgc_mappings.items()))
    json_data = [{
        "genome_ID": k,
        "BGC_ID": v
    } for k, v in genome_bgc_mappings.items()]
    json_data = {"mappings": json_data, "version": "1.0"}

    # validate json data
    validate(instance=json_data, schema=GENOME_BGC_MAPPINGS_SCHEMA)

    if output_file is None:
        output_file = bgc_dir / GENOME_BGC_MAPPINGS_FILENAME
    with open(output_file, "w") as f:
        json.dump(json_data, f)
    logger.info("Generated genome-BGC mappings file: %s", output_file)


def map_strain_to_bgc(strains: StrainCollection, bgcs: list[BGC]):
    """To set BGC object's strain with representative strain object.

    This method changes the list `bgcs` in place.

    Args:
        strains(StrainCollection): A collection of all strain objects.
        bgcs(list[BGC]): A list of BGC objects.

    Raises:
        KeyError: Strain id not found in the strain collection.
    """
    for bgc in bgcs:
        try:
            strain_list = strains.lookup(bgc.bgc_id)
            if len(strain_list) > 1:
                raise KeyError(
                    f"Multiple strain objects found for BGC id '{bgc.bgc_id}'."
                    f"BGC object accept only one strain."
                )
        except KeyError as e:
            raise KeyError(
                f"Strain id '{bgc.bgc_id}' from BGC object '{bgc.bgc_id}' "
                "not found in the strain collection.") from e
        bgc.strain = strain_list[0]


def map_bgc_to_gcf(bgcs: list[BGC], gcfs: list[GCF]):
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
                    "not found in the list of BGC objects.") from e
            gcf.add_bgc(bgc)


def filter_mibig_only_gcf(gcfs: list[GCF]) -> list[GCF]:
    """Filter out GCFs that contain only MIBiG BGC objects.

        This method returns a new list of GCFs that have at least one non-MIBiG
        BGC object as its child.
    """
    return [gcf for gcf in gcfs if gcf.has_mibig_only() is False]


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



@deprecated(version="1.3.3", reason="It is split to separate functions: " \
            "map_strain_to_bgc, map_bgc_to_gcf, filter_mibig_only_gcf, " \
            "get_bgcs_from_gcfs and get_strains_from_bgcs.")
def load_gcfs(bigscape_dir: str | PathLike, strains: StrainCollection,
              mibig_bgc_dict: dict[str, BGC], antismash_bgc_dict: dict[str,
                                                                       BGC],
              antismash_file_dict: dict[str, str], bigscape_cutoff: int):

    bigscape_dir = Path(bigscape_dir)
    product_class_cluster_file = bigscape_dir / "mix" / f"mix_clustering_c0.{bigscape_cutoff:02d}.tsv"
    network_annotations_file = bigscape_dir / "Network_Annotations_Full.tsv"

    new_bgc: BGC
    num_mibig: int = 0
    bgc_list: list[BGC] = []

    gcf_dict: dict[str, GCF] = {}
    gcf_list: list[GCF] = []

    used_strains: StrainCollection = StrainCollection()
    unknown_strains: dict[str, str] = {}

    # CG: bigscape data
    # parse the annotation files (<dataset>/bigscape/<cluster_name>/Network_Annotations_<cluster_name>.tsv
    # these contain fields:
    # - BGC name/ID [0]
    # - "Accession ID" [1]
    # - Description [2]
    # - Product prediction [3]
    # - Bigscape product type/class [4]
    # - Organism [5]
    # - Taxonomy [6]
    metadata = {}
    with open(network_annotations_file) as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)  # skip headers
        for line in reader:
            metadata[line[0]] = line

    # CG: bigscape data
    # "cluster files" are the various <class>_clustering_c0.xx.tsv files
    # - BGC name
    # - cluster ID
    with open(product_class_cluster_file, "rt") as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)  # skip headers
        for line in reader:
            bgc_name = line[0]
            family_id = line[1]

            # TODO: is it necessary to keep bigscape_class for GCF class?
            # get bgc annotations from bigscape file
            metadata_line = metadata[bgc_name]
            bigscape_class = metadata_line[4]

            # check strain
            try:
                strain = strains.lookup(bgc_name)
            except KeyError:
                logger.warning(f"Unknown strain ID: {bgc_name}")
                unknown_strains[bgc_name] = antismash_file_dict[bgc_name]
                continue

            # build new bgc
            if strain.id.startswith('BGC'):
                try:
                    new_bgc = mibig_bgc_dict[strain.id]
                except KeyError:
                    raise KeyError(f'Unknown MiBIG: {strain.id}')
                num_mibig += 1
            else:
                try:
                    new_bgc = antismash_bgc_dict[bgc_name]
                except KeyError:
                    raise KeyError(f'Unknown AntiSMASH BGC: {bgc_name}')

            new_bgc.strain = strain
            bgc_list.append(new_bgc)

            # build new gcf
            if family_id not in gcf_dict:
                new_gcf = GCF(family_id)
                gcf_dict[family_id] = new_gcf
                gcf_list.append(new_gcf)

            # link bgc to gcf
            gcf_dict[family_id].add_bgc(new_bgc)

            # add strain to used strains
            used_strains.add(strain)

    logger.info(
        '# MiBIG BGCs = {}, non-MiBIG BGCS = {}, total bgcs = {}, GCFs = {}, strains={}'
        .format(num_mibig,
                len(bgc_list) - num_mibig, len(bgc_list), len(gcf_dict),
                len(strains)))

    # filter out MiBIG-only GCFs)
    gcf_list, bgc_list, used_strains = _filter_gcfs(gcf_list, bgc_list,
                                                    used_strains)
    logger.info(
        '# after filtering, total bgcs = {}, GCFs = {}, strains={}, unknown_strains={}'
        .format(len(bgc_list), len(gcf_list), len(used_strains),
                len(unknown_strains)))

    return gcf_list, bgc_list, used_strains, unknown_strains


@deprecated(version="1.3.3", reason="It is split to separate functions: " \
            "filter_mibig_only_gcf, get_bgcs_from_gcfs and get_strains_from_bgcs.")
def _filter_gcfs(
    gcfs: list[GCF], bgcs: list[BGC], strains: StrainCollection
) -> tuple[list[GCF], list[BGC], StrainCollection]:
    """Remove a GCF from given GCF list if it only has MIBiG BGC members,
        correspondingly remove relevant BGC and strain from given list/collection.

        GCF and BGC internal id is updated to keep ids consectutive in a list.

    Args:
        gcfs(list[GCF]): list of GCF objects
        bgcs(list[BGC]): list of BGC objects
        strains(StrainCollection): StrainCollection object

    Returns:
        tuple[list[GCF], list[BGC], StrainCollection]: updated list of GCF
        objects, updated list of BGC objects and updated StrainCollection
        object.
    """
    gcfs_to_remove = set()
    bgcs_to_remove = set()

    for gcf in gcfs:
        num_non_mibig_bgcs = len(
            list(filter(lambda bgc: not bgc.is_mibig(), gcf.bgcs)))
        if num_non_mibig_bgcs == 0:
            gcfs_to_remove.add(gcf)
            for bgc in gcf.bgcs:
                bgcs_to_remove.add(bgc)

    for bgc in bgcs:
        if len(bgc.parents) == 0:
            bgcs_to_remove.add(bgc)

    for gcf in gcfs_to_remove:
        gcfs.remove(gcf)

    for bgc in bgcs_to_remove:
        bgcs.remove(bgc)
        if bgc.strain is not None:
            strains.remove(bgc.strain)

    logger.info(
        'Remove GCFs that has only MIBiG BGCs: removing {} GCFs and {} BGCs'.
        format(len(gcfs_to_remove), len(bgcs_to_remove)))

    return gcfs, bgcs, strains
