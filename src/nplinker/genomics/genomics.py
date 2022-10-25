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
import glob
import json
import os
import re
from Bio import SeqIO
from nplinker.logconfig import LogConfig
from nplinker.pairedomics import downloader
from nplinker.strains import Strain

from .bgc import BGC
from .gcf import GCF
from .mibigbgc import MiBIGBGC


logger = LogConfig.getLogger(__file__)

CLUSTER_REGION_REGEX = re.compile('(.+?)\\.(cluster|region)(\\d+).gbk$')


def parse_gbk_header(bgc):
    records = list(SeqIO.parse(bgc.antismash_file, format='gb'))
    if len(records) > 0:
        bgc.antismash_accession = records[0].name
        bgc.antismash_id = records[0].id


def loadBGC_from_cluster_files(strains, cluster_file_dict, ann_file_dict,
                               network_file_dict, mibig_bgc_dict,
                               mibig_json_dir, antismash_dir,
                               antismash_filenames, antismash_format,
                               antismash_delimiters):
    gcf_dict = {}
    gcf_list = []
    metadata = {}

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
    for a in ann_file_dict.values():
        with open(a) as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader)  # skip headers
            for line in reader:
                metadata[line[0]] = line

    num_mibig = 0
    num_missing_antismash = 0

    bgc_lookup = {}
    internal_bgc_id = len(mibig_bgc_dict)  # start numbering BGCs from here
    internal_gcf_id = 0

    bgc_list = [v for v in mibig_bgc_dict.values()]

    unknown_strains = {}

    logger.info(
        f'Using antiSMASH filename delimiters {antismash_delimiters}')

    # CG: bigscape data
    # "cluster files" are the various <class>_clustering_c0.xx.tsv files
    # - BGC name
    # - cluster ID
    for product_type, filename in cluster_file_dict.items():
        product_type = os.path.split(filename)[-1]
        product_type = product_type[:product_type.index('_')]
        with open(filename) as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader)  # skip headers
            for line in reader:
                name = line[0]
                family_id = int(line[1])
                if name.startswith('BGC'):
                    # removing the .<digit> suffix
                    nname = name[:name.index('.')]
                    strain = strains.lookup(nname)
                    if strain is None:
                        # if this happens, it probably means we have an MiBIG BGC which has been mistakenly
                        # excluded from the JSON database archive that NPLinker downloads. For more info
                        # see https://github.com/sdrogers/nplinker/issues/60#issuecomment-1086722952.
                        #
                        # To attempt to fix this issue without user intervention, try to download the
                        # missing BGC JSON data from the MiBIG website
                        if not downloader.download_mibig_bgc_json(
                                mibig_json_dir, nname):
                            # download failed, bail out here
                            raise Exception(
                                'Unknown MiBIG BGC: original={} / parsed={}'.
                                format(name, nname))
                        else:
                            # retrieved the file successfully but now have to parse it and add
                            # a new BGC to the existing set
                            strains, mibig_bgc_dict = append_mibig_library_json(
                                strains, mibig_bgc_dict, mibig_json_dir, nname,
                                internal_bgc_id)
                            logger.info(
                                'Appended MiBIG BGC {} to existing set'.format(
                                    nname))

                            # now can try again to lookup the strain, which should succeed this time
                            strain = strains.lookup(nname)
                            if strain is None:
                                # something is still wrong if this happens
                                raise Exception(
                                    'Unknown MiBIG BGC: original={} / parsed={}'
                                    .format(name, nname))

                            logger.info(
                                'MiBIG missing BGC workaround was successful')

                else:
                    parsednames = [
                        name[:name.index(d)] for d in antismash_delimiters
                        if name.find(d) != -1
                    ]
                    found = False
                    for parsedname in parsednames:
                        strain = strains.lookup(parsedname)
                        if strain is not None:
                            found = True
                            break
                        else:
                            # TODO hack to get crusemann working, should really update strain mappings?
                            for i in range(1, 3, 1):
                                tmp = f'{parsedname}.{i}'
                                strain = strains.lookup(tmp)
                                if strain is not None:
                                    found = True
                                    break
                            if found:
                                break

                    if not found:
                        logger.warning(
                            'Unknown strain ID: {} (from file {})'.format(
                                name, filename))
                        unknown_strains[name] = filename
                        continue

                # logger.debug('"{}" matched to {}'.format(name, strain))
                metadata_line = metadata[name]
                description = metadata_line[2]
                bigscape_class = metadata_line[4]
                product_prediction = metadata_line[3]

                # make a BGC object, reusing existing objects if they represent the same physical thing
                if not strain.id.startswith('BGC'):
                    if name in bgc_lookup:
                        new_bgc = bgc_lookup.get(name)
                    else:
                        # create a new BGC, increment internal ID and add to the list
                        new_bgc = BGC(internal_bgc_id, strain, name,
                                      product_prediction, description)
                        internal_bgc_id += 1
                        bgc_list.append(new_bgc)

                    if antismash_dir:
                        # TODO remove this at some point
                        if antismash_format == 'flat':
                            antismash_filename = os.path.join(
                                antismash_dir, new_bgc.name + '.gbk')
                            if not os.path.exists(antismash_filename):
                                logger.warn(
                                    '!!! Missing antismash file: {}'.format(
                                        antismash_filename))
                                num_missing_antismash += 1
                                # return None, None, None
                            new_bgc.set_filename(antismash_filename)
                            parse_gbk_header(new_bgc)
                        else:
                            new_bgc.set_filename(
                                antismash_filenames.get(new_bgc.name, None))
                            if new_bgc.antismash_file is None:
                                # TODO in some instances (e.g. with MSV000078836/Crusemann dataset),
                                # BiG-SCAPE output files appear to switch "cluster" with "region" in
                                # the lists of IDs. this leads to .gbk files not being matched up with
                                # BGCs even though they exist. Is this a sensible workaround???
                                if 'cluster' in new_bgc.name:
                                    new_bgc.set_filename(
                                        antismash_filenames.get(
                                            new_bgc.name.replace(
                                                'cluster', 'region'), None))
                                elif 'region' in new_bgc.name:
                                    new_bgc.set_filename(
                                        antismash_filenames.get(
                                            new_bgc.name.replace(
                                                'region', 'cluster'), None))

                            if new_bgc.antismash_file is None:
                                logger.warning(
                                    'Failed to find an antiSMASH file for {} {}'
                                    .format(new_bgc.name, new_bgc))
                                num_missing_antismash += 1
                            else:
                                parse_gbk_header(new_bgc)

                else:
                    num_mibig += 1
                    # TODO any reason not to supply the metadata fields that aren't set by
                    # make_mibig_bgc_dict since metadata_line is available here?
                    try:
                        new_bgc = mibig_bgc_dict[strain.id]
                    except KeyError:
                        raise Exception(f'Unknown MiBIG: {strain.id}')

                    new_bgc.description = description

                if family_id not in gcf_dict:
                    new_gcf = GCF(internal_gcf_id, family_id, bigscape_class)
                    gcf_dict[family_id] = new_gcf
                    gcf_list.append(new_gcf)
                    internal_gcf_id += 1

                gcf_dict[family_id].add_bgc(new_bgc)

                bgc_lookup[new_bgc.name] = new_bgc

    if num_missing_antismash > 0:
        logger.warn('{}/{} antiSMASH files could not be found!'.format(
            num_missing_antismash, len(bgc_list)))
        # print(list(antismash_filenames.keys()))

    logger.info(
        '# MiBIG BGCs = {}, non-MiBIG BGCS = {}, total bgcs = {}, GCFs = {}, strains={}'
        .format(num_mibig,
                len(bgc_list) - num_mibig, len(bgc_list), len(gcf_dict),
                len(strains)))

    # filter out irrelevant MiBIG BGCs (and MiBIG-only GCFs)
    bgc_list, gcf_list, strains = filter_mibig_bgcs(bgc_list, gcf_list,
                                                    strains)
    # update lookup table as well
    bgc_lookup = {bgc.name: bgc for bgc in bgc_list}

    logger.info(
        '# after filtering, total bgcs = {}, GCFs = {}, strains={}, unknown_strains={}'
        .format(len(bgc_list), len(gcf_list), len(strains),
                len(unknown_strains)))

    # load edge info - note that this should be done AFTER the filtering step above
    # so that it won't leave us with edges for BGCs that are no longer present
    logger.debug('Loading .network files')
    for filename in network_file_dict.values():
        with open(filename) as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader)  # skip headers
            # try to look up bgc IDs
            for line in reader:
                for i in range(2):
                    if line[i].startswith('BGC'):
                        # removing the .<digit> suffix
                        line[i] = line[i][:line[i].index('.')]

                if line[0] not in bgc_lookup or line[1] not in bgc_lookup:
                    # should indicate that one or both of these BGCs have been filtered out above
                    continue

                bgc_src = bgc_lookup[line[0]]
                bgc_dst = bgc_lookup[line[1]]
                bgc_src.edges.add(bgc_dst.id)

    return gcf_list, bgc_list, strains, unknown_strains


def filter_mibig_bgcs(bgcs, gcfs, strains):
    # remove the following MiBIG BGCs:
    # - parent set is empty (indicating never added to a GCF)
    # - any instances in a GCF with no other non-MiBIG BGCs
    to_remove_gcfs = set()
    to_remove_bgcs = set()
    for gcf in gcfs:
        # can ignore any with no MiBIGs
        if gcf.num_mibig_bgcs == 0:
            continue
        # now know this GCF has >0 MiBIG BGCs. Next step is to
        # check number of non-MiBIG BGCs, and if this is not at
        # least 1, throw away both GCF and BGC(s)
        if gcf.num_non_mibig_bgcs == 0:
            to_remove_gcfs.add(gcf)
            for bgc in gcf.bgcs:
                to_remove_bgcs.add(bgc)
                strains.remove(bgc.strain)
                bgc.parents.remove(gcf)

    for bgc in bgcs:
        if len(bgc.parents) == 0:
            strains.remove(bgc.strain)

    logger.info('Filtering MiBIG BGCs: removing {} GCFs and {} BGCs'.format(
        len(to_remove_gcfs), len(to_remove_bgcs)))

    # for GCFs just remove those that appear in to_remove_gcfs
    new_gcf_list = [gcf for gcf in gcfs if gcf not in to_remove_gcfs]
    # for BGCs do similar but also get rid of the objects never added to a GCF in the first place
    new_bgc_list = [
        bgc for bgc in bgcs
        if bgc not in to_remove_bgcs and len(bgc.parents) != 0
    ]

    # keep internal IDs consecutive
    for i in range(len(new_bgc_list)):
        new_bgc_list[i].id = i
        if len(new_bgc_list[i].parents) == 0:
            raise Exception(new_bgc_list[i])

    for i in range(len(new_gcf_list)):
        new_gcf_list[i].id = i

    return new_bgc_list, new_gcf_list, strains


def load_mibig_map(filename='mibig_gnps_links_q3_loose.csv'):
    mibig_map = {}
    with open(filename) as f:
        reader = csv.reader(f)
        next(reader)  # skip headers

        for line in reader:
            bgc = line[0]
            if bgc in mibig_map:
                mibig_map[bgc].append(line[3])
            else:
                mibig_map[bgc] = [line[3]]
    return mibig_map


def load_mibig_library_json(mibig_json_directory):
    mibig = {}
    files = glob.glob(mibig_json_directory + os.sep + '*.json')
    logger.info(f"Found {len(files)} MiBIG json files")
    for file in files:
        with open(file) as f:
            bgc_id = file.split(os.sep)[-1].split('.')[0]
            mibig[bgc_id] = json.load(f)
    return mibig


def extract_mibig_json_data(data):
    if 'general_params' in data:
        accession = data['general_params']['mibig_accession']
        biosyn_class = data['general_params']['biosyn_class'][0]
    else:  # 2.0(+)
        accession = data['cluster']['mibig_accession']
        biosyn_class = data['cluster']['biosyn_class'][0]

    return accession, biosyn_class


def append_mibig_library_json(strains, mibig_bgc_dict, mibig_json_directory,
                              bgc_id, internal_id):
    json_data = json.load(
        open(os.path.join(mibig_json_directory, f'{bgc_id}.json'),
             'rb'))
    accession, biosyn_class = extract_mibig_json_data(json_data)
    strain = Strain(accession)
    new_bgc = MiBIGBGC(internal_id, strain, accession, biosyn_class)
    mibig_bgc_dict[accession] = new_bgc
    strains.add(strain)
    return strains, mibig_bgc_dict


def make_mibig_bgc_dict(strains, mibig_json_directory, version):
    mibig_dict = load_mibig_library_json(mibig_json_directory)
    mibig_bgc_dict = {}
    i = 0
    for name, data in list(mibig_dict.items()):
        accession, biosyn_class = extract_mibig_json_data(data)
        strain = Strain(accession)
        new_bgc = MiBIGBGC(i, strain, accession, biosyn_class)
        mibig_bgc_dict[accession] = new_bgc
        strains.add(strain)
        i += 1
    return mibig_bgc_dict
