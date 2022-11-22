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

import csv
import os
import re
from Bio import SeqIO
from nplinker.genomics.mibig import MibigBGC
from nplinker.logconfig import LogConfig
from nplinker.strain_collection import StrainCollection
from .bgc import BGC
from .gcf import GCF

logger = LogConfig.getLogger(__name__)

CLUSTER_REGION_REGEX = re.compile('(.+?)\\.(cluster|region)(\\d+).gbk$')


def parse_gbk_header(bgc):
    """Read AntiSMASH BGC .gbk file to get BGC name and id"""
    records = list(SeqIO.parse(bgc.antismash_file, format='gb'))
    if len(records) > 0:
        bgc.antismash_accession = records[0].name
        bgc.antismash_id = records[0].id


def load_gcfs(strains: StrainCollection, product_class_cluster_file: str,
              network_annotations_file: str,
              mibig_bgc_dict: dict[str, MibigBGC],
              antismash_bgc_dict: dict[str, BGC],
              antismash_file_dict: dict[str, str]):
    metadata = {}

    num_mibig = 0
    internal_bgc_id = 0
    bgc_list = []

    internal_gcf_id = 0
    gcf_dict = {}
    gcf_list = []

    used_strains = StrainCollection()
    unknown_strains = {}

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
            family_id = int(line[1])

            # get bgc annotations from bigscape file
            metadata_line = metadata[bgc_name]
            description = metadata_line[2]
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
                new_bgc.description = description
            else:
                try:
                    new_bgc = antismash_bgc_dict[bgc_name]
                except KeyError:
                    raise KeyError(f'Unknown AntiSMASH BGC: {bgc_name}')
                # update strain to make sure consistent strain id
                new_bgc.strain = strain

            new_bgc.id = internal_bgc_id
            bgc_list.append(new_bgc)
            internal_bgc_id += 1

            # build new gcf
            if family_id not in gcf_dict:
                new_gcf = GCF(internal_gcf_id, family_id, bigscape_class)
                gcf_dict[family_id] = new_gcf
                gcf_list.append(new_gcf)
                internal_gcf_id += 1

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
    gcf_list, bgc_list, used_strains = _filter_gcfs(gcf_list, bgc_list, used_strains)
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
        if gcf.num_non_mibig_bgcs == 0:
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

    # keep internal IDs consecutive
    for index, bgc in enumerate(bgcs):
        bgc.id = index
    for index, gcf in enumerate(gcfs):
        gcf.id = index

    logger.info(
        'Remove GCFs that has only MIBiG BGCs: removing {} GCFs and {} BGCs'.
        format(len(gcfs_to_remove), len(bgcs_to_remove)))

    return gcfs, bgcs, strains
