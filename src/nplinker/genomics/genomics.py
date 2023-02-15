# Copyright 2021 The NPLinker Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
from __future__ import annotations
import csv
import re
from Bio import SeqIO
from nplinker.logconfig import LogConfig
from nplinker.strain_collection import StrainCollection
from .bgc import BGC
from .gcf import GCF
from os import PathLike
import os
from pathlib import Path

logger = LogConfig.getLogger(__name__)

CLUSTER_REGION_REGEX = re.compile('(.+?)\\.(cluster|region)(\\d+).gbk$')


def load_gcfs(bigscape_dir: str | PathLike,
              strains: StrainCollection,
              mibig_bgc_dict: dict[str, BGC],
              antismash_bgc_dict: dict[str, BGC],
              antismash_file_dict: dict[str, str],
              bigscape_cutoff: int):

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
                new_gcf = GCF(family_id, bigscape_class)
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
        num_non_mibig_bgcs = len(list(filter(lambda bgc: not bgc.is_mibig(), gcf.bgcs)))
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
        strains.remove(bgc.strain)

    logger.info(
        'Remove GCFs that has only MIBiG BGCs: removing {} GCFs and {} BGCs'.
        format(len(gcfs_to_remove), len(bgcs_to_remove)))

    return gcfs, bgcs, strains
