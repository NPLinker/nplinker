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
import json
import logging
import os
import re


logger = logging.getLogger(__name__)

# this will match strings like '...cluster001.gbk' or '...region022.gbk',
# and allow the number to be extracted easily
CLUSTER_REGION_REGEX = re.compile("(.+?)\\.(cluster|region)(\\d+).gbk$")

# This module contains a couple of parsers for antiSMASH knownclusterblast output
# files: one for the newer per-directory .JSON blob and one for the legacy
# format with one text file per gbk.


class KCBJSONParser:
    """Parser for the large .json files antiSMASH generates as part of its output.

    This is supposed to do the same job as the KCBTextParser class without relying
    on parsing the legacy-format text files.

    Args:
        bgcs: a list of NPLinker BGC objects created from the same antiSMASH output folder
    """

    def __init__(self, bgcs):
        # check all the linked files exist
        for bgc in bgcs:
            if not os.path.exists(bgc.antismash_file):
                raise Exception('KCBJSONParser failed to find file "{}"'.format(bgc.antismash_file))

        logger.info(f"KCBJSONParser({len(bgcs)} BGCs)")

        # find the JSON file: TODO is the assumption of there only being a single .json
        # file always going to work? otherwise have to try guessing the name based on
        # genome IDs
        prefix = os.path.dirname(bgcs[0].antismash_file)
        json_files = list(filter(lambda f: f.endswith(".json"), os.listdir(prefix)))
        logger.info("Found {} JSON files in {}".format(len(json_files), prefix))

        if len(json_files) == 0:
            logger.warning("Unable to find an antiSMASH JSON output file in {}".format(prefix))
            self.json_filename = None
            return

        self.json_filename = os.path.join(prefix, json_files[0])
        logger.info(f"Using JSON file {self.json_filename}")

    def parse_hits(self):
        if self.json_filename is None:
            return None

        self.collected_hits = {}

        # NOTE: since some of these JSON files are ridiculously large (>50MB in some cases)
        # it can cause OOM errors (with the stdlib json module and a couple of other
        # 3rd party ones). Best way to avoid this is by using 64-bit Python.
        # TODO: catch exception and print message recommending 64-bit
        logger.info("Loading antiSMASH JSON data from {}".format(self.json_filename))
        data = json.load(open(self.json_filename))

        # JSON data structure: depending on the gbks in the source folder, the
        # structure of the results may be different from one to the next. from
        # looking at some examples, what should be expected is a structure like
        # this outline ([] = array/list object):
        #   <top>
        #       [records]
        #           id
        #           name
        #           description
        #           [features]
        #           [modules]
        #               ...
        #               [antismash.modules.clusterblast] (may be missing)
        #                   record_id
        #                   [knowncluster]
        #                       record_id
        #                       [results]
        #                           region_number
        #                           total_hits
        #                           [ranking] (may be empty if total_hits = 0)
        #                               (tuple_pt1, tuple_pt2)
        #                               ...
        #                       mibig_entries
        #                           (lists of IDs named by region numbers)
        #
        # every gbk with a distinct prefix maps to an entry in the [records] list.
        # however there are also entries for gbks which don't exist in the source
        # folder, don't know why. These can easily be excluded by checking for
        # the presence of the "antismash.modules.clusterblast" module results in
        # the [modules] list for that record.
        #
        # the next step is to access the [results] list nested further down in
        # the structure. the number of elements in this list are dependent on the
        # number of gbks that are different regions of the same genome. e.g. if
        # you have foobar.region001.gbk and foobar.region002.gbk, there should be
        # a single 'record' for 'foobar' and a 'result' for each region. The
        # region_number field can easily be used to distinguish these.
        #
        # finally there is the 'ranking' list. this is where all the data we need
        # is located. it has a weird structure where every entry is a 2-tuple with
        # this content:
        #   1st element: has fields accession, cluster_label, description, and
        #       cluster_type, plus lists 'proteins' and 'tags'
        #   2nd element: has fields hits, blast_score, core_gene_hits,
        #       synteny_score, core_bonus, and a 'pairings' list
        #
        # each element of the 'pairings' list is in turn a 3-tuple:
        #   1st element: a |-delimited field with multiple values in it. a typical
        #       example is 'input|c1|0-2665|-|B098_RS0100005|indolepyruvate'
        #   2nd element: a single integer, which seems to be the 0-based index
        #       of this 'pairing' in the list?
        #   3rd element: an object with fields name, genecluster, start, end,
        #       strand, annotation, perc_ident, blastscore, perc_coverage, evalue,
        #       and locus_tag
        #
        #  to correctly parse this, need to iterate over the list of records, and
        #  for each one that has clusterblast results iterate over those, matching
        #  things up to the BGC objects created by nplinker as we go along.

        for i, rec in enumerate(data["records"]):
            hits = self._parse(rec)
            if hits is not None:
                self.collected_hits.update(hits)

        logger.info(
            "KCBJSONParser: collected {} total hit entries".format(len(self.collected_hits))
        )

        return self.collected_hits

    def _parse(self, record):
        """Parses the knownclusterblast data for a single 'record' entry."""
        modules = record.get("modules", None)
        if modules is None or "antismash.modules.clusterblast" not in modules:
            # this probably isn't an error, the JSON often seems to contain entries
            # for gbks that don't exist in the source folder (??)
            return None

        kcb = modules["antismash.modules.clusterblast"]["knowncluster"]
        record_id = kcb["record_id"]

        # each 'record' may contain multiple results, and each 'result' may be
        # linked to a different BGC and region number, so need to keep track of
        # all this to match them to NPLinker BGC objects later
        record_hits = {}

        for i, result in enumerate(kcb["results"]):
            # check number of rankings, can ignore if zero
            if len(result["ranking"]) == 0:
                continue

            # step 1: extract the keys of the 'mibig_entries' dict for the current region number (there should be
            # the same number of 'results' as 'mibig_entries'). There appear to be fewer entries in these lists than
            # appear in the text files where both exist, but I can't find another source for them and the extras
            # don't appear to matter much as they aren't referenced later in the text files or by the original parser
            gene_ids = list(kcb["mibig_entries"][str(i + 1)].keys())

            # step 2: extract the MiBIG BGC IDs + product names (=parsing the Significant hits table in text file)
            sig_hits = [
                (entry[0]["accession"], entry[0]["description"]) for entry in result["ranking"]
            ]

            # give up if there aren't any of these
            if len(sig_hits) == 0:
                continue

            # step 3: extract BGC IDs (=parsing of "... subject cluster" tables in text file)
            # (this is a list of lists)
            all_mibig_genes = [entry[0]["proteins"] for entry in result["ranking"]]

            region_number = int(result["region_number"])
            # should never have two identical region numbers for the same ID
            assert region_number not in record_hits

            # want to construct the same data structure as the original parser
            # should contain:
            #   - all_bgc_genes => gene_ids
            #   - all_mibig_genes => all_mibig_genes
            #   - invididual_hits => list populated based on blast_hits
            hit = {}
            assert len(all_mibig_genes) == len(result["ranking"])

            # step 4: extract the blast hit information from the 'rankings' and
            # 'pairings' lists (=parsing of "Table of Blast hits" in text file)

            # each "hit" shares these
            hit = {"all_bgc_genes": gene_ids}
            individual_hits = []
            for j, ranking in enumerate(result["ranking"]):
                hit["all_mibig_genes"] = all_mibig_genes[j]
                hit["mibig_id"] = sig_hits[j][0]
                for pairing in ranking[1]["pairings"]:
                    tc1 = pairing[0].split("|")[4]
                    individual_hits.append(
                        {
                            "source_bgc_gene": tc1,
                            "mibig_bgc_gene": pairing[2]["name"],
                            "identity_percent": int(pairing[2]["perc_ident"]),
                            "blast_score": int(pairing[2]["blastscore"]),
                        }
                    )

            hit["individual_hits"] = individual_hits
            record_hits[region_number] = hit
            # logger.debug('\tSig hit {}/{} for {}: found {} individual hits'.format(j+1, len(sig_hits), sig_hits[j][0], len(individual_hits)))

        # don't create entries with no results
        if len(record_hits) == 0:
            return {}

        return {record_id: record_hits}


class KCBTextParser:
    """Parser for antismash knownclusterblast text output files."""

    def __init__(self, filename):
        if not os.path.exists(filename):
            raise Exception(f'KCBTextParser failed to find file "{filename}"')

        self.bgc_genes = set()
        self.mibig_bgcs = []
        self.hits = {}

        # Structure outline for quick reference:
        # <top of file>
        # Clusterblast scores for <ID>
        #
        # Table of genes, locations, strands and annotations of query cluster:
        # <multiple lines in table> (this is top_block below)
        #
        # Significant hits:
        # <multiple hits> (this is second_block below)
        #
        # Details:
        #
        # >>
        # <multiple lines of details for hit 1>
        #
        # >>
        # <multiple lines of details for hit 2>
        #
        # (continues up to hit n)
        # EOF

        with open(filename) as f:
            line = next(f)
            while not line.startswith("Table of genes"):
                line = next(f)
            # now we're in the first block. this is the section beginning with the
            # line "Table of genes, locations, strands and annotations of query cluster:"
            # and ending at "Significant hits:"
            top_block = []
            while True:
                line = next(f)
                if line.startswith("Significant"):
                    break
                else:
                    if len(line) > 1:
                        top_block.append(line.rstrip())

            # now we're in the second block. this is the section beginning with "Significant hits:"
            # and ending at "Details:"
            second_block = []
            while True:
                line = next(f)
                if line.startswith("Details"):
                    break
                else:
                    if len(line) > 1:
                        second_block.append(line.rstrip())

            # find the start of the first ">>" delimited section
            while True:
                try:
                    line = next(f)
                    if line.startswith(">>"):
                        break
                except Exception:
                    # EOF
                    return

            details = []
            finished = False
            while not finished:
                temp_list = []
                while True:
                    try:
                        line = next(f)
                        if line.startswith(">>"):
                            # finished current section
                            details.append(temp_list)
                            break
                        else:
                            if len(line) > 1:
                                temp_list.append(line.rstrip())
                    except Exception:
                        details.append(temp_list)
                        finished = True
                        break

            # do some processing on the blocks
            # firstly, extract the genes from the BGC -- stored in the first block
            # these are the gene names in the BGC
            # e.g. given a line "STROP_RS12480	2788559	2789786	-	cytochrome P450	"
            # this would add "STROP_RS12480" to self.bgc_genes
            for line in top_block:
                tokens = line.split()
                self.bgc_genes.add(tokens[0])

            # secondly, extract the MiBIG BGCs that are mentioned here
            # e.g. given a line "1. BGC0000279_c1	Xantholipin_biosynthetic_gene_cluster"
            # this would take "BGC0000279_c1" as the bgc_id and "Xantholipin_biosynthetic_gene_cluster"
            # as the bgc_product_name
            for line in second_block:
                tokens = line.split()
                bgc_id = tokens[1]
                bgc_product_name = tokens[2]
                self.mibig_bgcs.append((bgc_id, bgc_product_name))

            # details is a list of lists, where each sublist contains the set of lines
            # for one ">>" delimited section giving the detailed information on a
            # particular hit
            for i, detail in enumerate(details):
                # the first line is of the form "1. BGC0000279_c1"
                current_bgc_id = detail[0].split()[1]  # this is the MiBIG ID, e.g. BGC0001666

                self.hits[current_bgc_id] = {}
                self.hits[current_bgc_id]["all_bgc_genes"] = self.bgc_genes

                # if the detail sections don't appear in the same order as the "Significant hits"
                # section then something is wrong and can't continue parsing this
                if current_bgc_id != self.mibig_bgcs[i][0]:
                    raise Exception(
                        "Mismatched ordering found in knownclusterblast file {}".format(filename)
                    )

                # get the line numbers where these two tables start
                table_pos = detail.index(
                    "Table of genes, locations, strands and annotations of subject cluster:"
                )
                pos = detail.index(
                    "Table of Blast hits (query gene, subject gene, %identity, blast score, %coverage, e-value):"
                )

                # get the list of all of the genes within the MiBIG BGC by iterating over
                # the lines in the "Table of genes ..." table
                all_mibig_genes = []
                for line in detail[table_pos + 1 : pos]:
                    all_mibig_genes.append(line.split()[0])

                if len(all_mibig_genes) == 0:
                    logger.warning(
                        "KCBTextParser failed to extract any MiBIG genes from file {}, BGC ID {}".format(
                            filename, current_bgc_id
                        )
                    )
                    # just continue to the next section here because this currently makes the rest of the processing invalid
                    del self.hits[current_bgc_id]
                    continue

                self.hits[current_bgc_id]["all_mibig_genes"] = all_mibig_genes
                self.hits[current_bgc_id]["individual_hits"] = []

                for line in detail[pos + 1 :]:
                    tokens = line.split()
                    bgc_id = tokens[0]
                    self.hits[current_bgc_id]["individual_hits"].append(
                        {
                            "source_bgc_gene": tokens[0],
                            "mibig_bgc_gene": tokens[1],
                            "identity_percent": int(tokens[2]),
                            "blast_score": int(tokens[3]),
                        }
                    )

    @staticmethod
    def get_kcb_filename_from_bgc(bgc):
        """Given a BGC object, return the filename of the corresponding knownclusterblast .txt file (if any).

        This method attempts to derive the name of the knownclusterblast output file for the
        supplied BGC object, using the original path + filename of the .gbk that the BGC was
        sourced from during the loading process.

        Assumptions made:
            - the .gbk files are grouped in subdirectories with <dataset>/antismash/
            - each of these subdirectories contains a "knownclusterblast" subdirectory
            - within that, there is a single .txt file for each .gbk

        For some reason there doesn't appear to be a simple way of matching up the .gbk/.txt
        files using their contents, only the filenames. And there are various different
        possibilities for the naming schemes, so this method tries to account for the
        known variants discovered so far.

        Args:
            bgc: a BGC object

        Returns:
            A string containing the filename of the knownclusterblast .txt file for the
            supplied BGC, or None if an error occurred/file doesn't exist
        """
        if bgc.antismash_file is None:
            logger.warning(f"BGC {bgc} has no antismash_file set")
            return None

        # expecting to find the .txt files inside a 'knownclusterblast' subdir in the
        # same location as the .gbk file, give up if that doesn't exist
        base_path = os.path.join(os.path.dirname(bgc.antismash_file), "knownclusterblast")

        if not os.path.exists(base_path):
            logger.warning(
                'Expected "knownclusterblast" directory not found at "{}"'.format(base_path)
            )
            return None

        # get the name of the .gbk file itself (no path)
        genbank_file = os.path.split(bgc.antismash_file)[1]

        if "region" in genbank_file:
            # case 1: assume the genbank files have a <someID>.region<num>.gbk naming scheme.
            # this seems to be the case for all the more recent datasets i've seen.

            # use a regex to extract the portions of the filename we're interested in
            regex_obj = CLUSTER_REGION_REGEX.search(genbank_file)
            number = int(regex_obj.group(3))
            prefix = regex_obj.group(1)

            # construct the expected filename using the prefix and number
            # (note no leading zeroes on the number)
            kcb_name = os.path.join(base_path, f"{prefix}_c{number}.txt")
        elif "cluster" in genbank_file:
            # case 2: assume the genbank files have a <someID>.cluster<num>.gbk naming scheme.
            # this is the case with the Crusemann dataset on the paired platform among others.

            # use a regex to extract the number from the filename
            regex_obj = CLUSTER_REGION_REGEX.search(genbank_file)
            number = int(regex_obj.group(3))

            # construct the expected filename
            # (note no leading zeroes on the number)
            kcb_name = os.path.join(base_path, f"cluster{number}.txt")
        else:
            logger.warning(
                "Unknown GenBank file naming scheme, failed to determine knownclusterblast filenames!"
            )
            return None

        return kcb_name
