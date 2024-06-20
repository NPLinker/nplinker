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
import logging
import os
from nplinker.defaults import NPLINKER_APP_DATA_DIR
from nplinker.scoring.rosetta.rosetta_hit import RosettaHit
from ...genomics import BGC
from ...parsers.kcb import KCBJSONParser
from ...parsers.kcb import KCBTextParser
from ...pickler import load_pickled_data
from ...pickler import save_pickled_data
from .spec_lib import SpecLib


logger = logging.getLogger(__name__)


class Rosetta:
    DEF_MS1_TOL = 100
    DEF_MS2_TOL = 0.2
    DEF_SCORE_THRESH = 0.5
    DEF_MIN_MATCH_PEAKS = 1

    PARAM_VERSION = 1

    def __init__(self, nplinker, ignore_genomic_cache=False):
        self._nplinker = nplinker
        self._mgf_data = {}
        self._csv_data = {}
        self._mgf_path = os.path.join(NPLINKER_APP_DATA_DIR, "matched_mibig_gnps_update.mgf")
        self._csv_path = os.path.join(NPLINKER_APP_DATA_DIR, "matched_mibig_gnps_update.csv")
        self._data_path = NPLINKER_APP_DATA_DIR
        self._root_path = nplinker.root_dir
        self._dataset_id = nplinker.dataset_id
        self._ignore_genomic_cache = ignore_genomic_cache
        self._pickle_dir = os.path.join(nplinker.root_dir, "rosetta")
        if not os.path.exists(self._pickle_dir):
            os.makedirs(self._pickle_dir, exist_ok=True)
        self._speclib_pickle_path = os.path.join(self._pickle_dir, "SpecLib.pckl")
        self._spechits_pickle_path = os.path.join(self._pickle_dir, "spec_hits.pckl")
        self._bgchits_pickle_path = os.path.join(self._pickle_dir, "bgc_hits.pckl")
        self._rhits_pickle_path = os.path.join(self._pickle_dir, "RosettaHits.pckl")
        self._params_pickle_path = os.path.join(self._pickle_dir, "RosettaParams.pckl")

        if not os.path.exists(self._mgf_path):
            logger.warning(
                "Failed to load Rosetta data ({}), matching disabled".format(self._mgf_path)
            )
            return
        if not os.path.exists(self._csv_path):
            logger.warning(
                "Failed to load Rosetta data ({}), matching disabled".format(self._csv_path)
            )
            return

        self._gnps2mibig = None
        self._mibig2gnps = None
        self._mibig2bgc = {}

        self.speclib = None
        self._rosetta_hits = []
        self._spec_hits = None
        self._bgc_hits = None

    @property
    def bgc_hits(self):
        return self._bgc_hits

    @property
    def spec_hits(self):
        return self._spec_hits

    def _load_csv(self, csv_path):
        logger.info("constructing rosetta dicts")

        self._gnps2mibig = {}
        self._mibig2gnps = {}

        with open(csv_path) as f:
            rdr = csv.reader(f, delimiter=",")
            headers = next(rdr)
            for line in rdr:
                gnps, mibig = line[0], line[3]
                if gnps in self._gnps2mibig:
                    self._gnps2mibig[gnps].append(mibig)
                else:
                    self._gnps2mibig[gnps] = [mibig]

                if mibig in self._mibig2gnps:
                    self._mibig2gnps[mibig].append(gnps)
                else:
                    self._mibig2gnps[mibig] = [gnps]

    def _generate_spec_hits(self, spectra, ms1_tol, ms2_tol, score_thresh, min_match_peaks):
        spec_hits = {}
        for i, sp in enumerate(spectra):
            hits = self.speclib.spectral_match(
                sp,
                ms2_tol=ms2_tol,
                min_match_peaks=min_match_peaks,
                ms1_tol=ms1_tol,
                score_thresh=score_thresh,
            )
            if len(hits) > 0:
                spec_hits[sp] = hits
            if i % 100 == 0:
                logger.info("Searching for spectral hits {}/{}".format(i, len(spectra)))

        save_pickled_data(spec_hits, self._spechits_pickle_path)
        return spec_hits

    def _generate_speclib(self):
        logger.warning("No pickle SpecLib found, generating (this will take some time!)...")
        self.speclib = SpecLib(self._mgf_path)
        self.speclib._load_mgf()
        self.speclib.filter()

        logger.info("Finished generating SpecLib")

        save_pickled_data(self.speclib, self._speclib_pickle_path)

    def _generate_bgc_hits(self, bgcs):
        self._bgc_hits = {}
        errors = 0

        # this method is a bit messy because it tries to handle a couple of different
        # routes to extracting knownclusterblast results:
        #
        #   - newer antiSMASH runs should produce a single .JSON file in each directory
        #       of .gbk files. this one file can be parsed to extract all the info that the
        #       rosetta code requires for all of the .gbks (uses KCBJSONParser)
        #   - older datasets will only include the now-deprecated text format results, which
        #       should be located in a subdirectory of each of the .gbk directories called
        #       "knownclusterblast". in these instances there will be a single .txt file
        #       for each .gbk
        #   - some datasets may have no knownclusterblast results available
        #
        # it's also possible for all 3 of these to turn up within the same dataset. for example
        # when downloading from the paired omics platform, the antiSMASH data is downloaded
        # separately for each genome from the antiSMASH DB and this means the format may not
        # be consistent.

        logger.debug("Collecting BGC hit information...")
        # go through the list of all available BGCs (ignoring MiBIGBGC instances) and
        # group them by the directory they appear in
        bgc_groups = {}
        skipped, errors = 0, 0

        for bgc in bgcs:
            if bgc.antismash_file is None or isinstance(bgc, BGC):
                skipped += 1
                continue

            prefix = os.path.dirname(bgc.antismash_file)
            if prefix not in bgc_groups:
                bgc_groups[prefix] = [bgc]
            else:
                bgc_groups[prefix].append(bgc)

        logger.debug("{} BGC groups based on filenames".format(len(bgc_groups)))

        for prefix, prefix_bgcs in bgc_groups.items():
            logger.debug(
                "Attempting to parse JSON data for prefix {} with {} BGCs".format(
                    prefix, len(prefix_bgcs)
                )
            )
            # preferred option is to parse the results for the whole group using the
            # JSON file, but this may not be available...
            json_hits = KCBJSONParser(prefix_bgcs).parse_hits()
            matched_bgcs = {}

            if json_hits is not None:
                # number of hits can often be less than number of BGCs (e.g. if no significant hits found)
                sum_hits = sum(len(json_hits[x]) for x in json_hits)
                logger.debug(
                    "JSON parsing was successful! Returned {} hits from {} BGCs".format(
                        sum_hits, len(prefix_bgcs)
                    )
                )

                # unlike the KCBTextParser where each set of results is easy to link
                # back to the appropriate BGC object, here we need to do some extra work
                # to ensure we have everything matched up correctly.
                #
                # seems like the best way to do this is to rely on the BGC attributes
                # parsed from the .gbks during the loading process (including region
                # numbers) as these should match up directly with the JSON data.

                for pbgc in prefix_bgcs:
                    # the "normal" case appears to be that you'll have a directory
                    # containing multiple gbks with the same accession and different
                    # region numbers, e.g. ABC123.region001, ABC123.region002, ...
                    # and should be able to expect that every "hit" comes from a
                    # .gbk that exists in the directory.
                    #
                    # in ideal circumstances, the json_hits structure will end up
                    # containing a single top level accession (ABC123) and then
                    # at the next level down one region number for each of the
                    # BGCs with signficant hits. this makes matching BGC objects
                    # quite simple.
                    #
                    # however in other cases there appear to be a mix of accession
                    # IDs in the same antiSMASH directory. so you can have collections
                    # where the filenames are e.g. ABC123.region001, DEF456.region001,
                    # GHI789.region001, ... (including multiple regions for the same
                    # accession). this is more difficult to handle because the JSON
                    # data isn't a direct match with that parsed from the .gbks
                    # themselves. for example, the gbk might report a region number of
                    # 7 while the corresponding JSON result has a region number of 1.
                    # possible workaround is to take the gbk region number and check
                    # if it appears in the filename of the gbk???

                    if pbgc.antismash_id in json_hits:
                        # simplest case where there's a direct match on region number
                        if pbgc.antismash_region in json_hits[pbgc.antismash_id]:
                            logger.debug(
                                "Matched {} using {} + region{:03d}!".format(
                                    pbgc.antismash_file, pbgc.antismash_id, pbgc.antismash_region
                                )
                            )
                            if pbgc not in matched_bgcs:
                                matched_bgcs[pbgc] = {}

                            hit = json_hits[pbgc.antismash_id][pbgc.antismash_region]
                            matched_bgcs[pbgc][hit["mibig_id"]] = hit
                            continue
                        else:
                            # if the above case doesn't apply, check through every
                            # region number available for the antismash ID we have,
                            # and check if the original filename contains that region
                            # number. if so assume it is the correct match.
                            for region in json_hits[pbgc.antismash_id]:
                                if pbgc.antismash_file.endswith(
                                    "region{:03d}.gbk".format(pbgc.antismash_region)
                                ):
                                    logger.debug(
                                        "Matched {} using fallback {} + region{:03d} (orig={})".format(
                                            pbgc.antismash_file,
                                            pbgc.antismash_id,
                                            region,
                                            pbgc.antismash_region,
                                        )
                                    )
                                    if pbgc not in matched_bgcs:
                                        matched_bgcs[pbgc] = {}

                                    hit = json_hits[pbgc.antismash_id][region]
                                    matched_bgcs[pbgc][hit["mibig_id"]] = hit
                                    break
                    else:
                        # this could simply mean no significant hits found
                        logger.info(
                            "Found no matching hits for BGC ID={}, region={}, file={}".format(
                                pbgc.antismash_id, pbgc.antismash_region, pbgc.antismash_file
                            )
                        )

            else:
                # ... if JSON parsing failed, fall back to the original text parser. this
                # must be called on each BGC individually
                logger.debug("JSON parsing failed, falling back to text instead")
                for i, bgc in enumerate(prefix_bgcs):
                    kcb_name = KCBTextParser.get_kcb_filename_from_bgc(bgc)
                    if kcb_name is not None:
                        try:
                            parser = KCBTextParser(kcb_name)
                            if len(parser.hits) > 0:
                                matched_bgcs[bgc] = parser.hits
                        except Exception as e:
                            logger.warning(e)
                            errors += 1

            logger.debug("Found matches for {}/{} bgcs".format(len(matched_bgcs), len(prefix_bgcs)))
            if len(matched_bgcs) != len(prefix_bgcs):
                # not necessarily fatal but probably not good either
                logger.warning(
                    "Failed to match {} BGCs to hits in directory {}!".format(
                        len(prefix_bgcs) - len(matched_bgcs), prefix
                    )
                )

            # now insert the matched hits into the _bgc_hits structure so the original code
            # below can parse them in the same way as the text parser results
            self._bgc_hits.update(matched_bgcs)

        # make reverse dict
        self._mibig2bgc = {}
        for bgc, hits in self._bgc_hits.items():
            for mibig_bgc_id in hits:
                if mibig_bgc_id not in self._mibig2bgc:
                    self._mibig2bgc[mibig_bgc_id] = set()
                self._mibig2bgc[mibig_bgc_id].add(bgc)

        logger.info(f"Completed, {len(self._bgc_hits)} BGC hits found")
        if errors > 0:
            logger.warning(
                "Some knownclusterblast files could not be loaded, results may be incomplete"
            )
        save_pickled_data((self._bgc_hits, self._mibig2bgc), self._bgchits_pickle_path)

    def _collect_rosetta_hits(self):
        self._rosetta_hits = []
        bgc_summary_scores = self.generate_bgc_summary_scores()
        for spec, data in self._spec_hits.items():
            for gnps_id, score in data:
                for mibig_id in self._gnps2mibig[gnps_id]:
                    if mibig_id in self._mibig2bgc:
                        for bgc in self._mibig2bgc[mibig_id]:
                            # get the bgc score
                            bgc_score = bgc_summary_scores[bgc][mibig_id]
                            self._rosetta_hits.append(
                                RosettaHit(spec, gnps_id, mibig_id, bgc, score, bgc_score)
                            )
        logger.info(f"Found {len(self._rosetta_hits)} rosetta hits!")

    def generate_bgc_summary_scores(self):
        # process the hit to compress it into more useful info
        # computes the total identity score for each mibig entry
        # and divides by the number of mibig genes
        # i.e. the score represents how much of the mibig is
        # reflected in the source bgc

        # we have one entry per mibig that is linked to this BGC and
        # a single score (that will vary between 0 and 1)

        processed = {}
        # bgc_hits = {BGC object: list of hits from KCB}
        # each list is a list of dicts
        for bgc, hit in self._bgc_hits.items():
            mibig_bgcs = list(hit.keys())
            scores = {}
            for mibig_id in mibig_bgcs:
                n_mibig_genes = len(hit[mibig_id]["all_mibig_genes"])
                if n_mibig_genes == 0:
                    logger.warning(
                        "Found a BGC entry with zero genes, this should never happen! (BGC={})".format(
                            bgc
                        )
                    )
                    continue
                # n_source_genes = len(hit[mibig_id]['all_bgc_genes'])
                total_hit_identity = 0
                for hit_gene in hit[mibig_id]["individual_hits"]:
                    identity_percent = hit_gene["identity_percent"]
                    total_hit_identity += identity_percent / 100.0
                score = total_hit_identity / n_mibig_genes
                scores[mibig_id] = score
            processed[bgc] = scores
        return processed

    def run(self, spectra, bgcs, ms1_tol, ms2_tol, score_thresh, min_match_peaks):
        """Function which actually computes the rosetta score somehow."""
        params_ok = self._load_cached_params(ms1_tol, ms2_tol, score_thresh, min_match_peaks)

        # if any parameters have been changed or version mismatch found, delete all cached files
        if not params_ok:
            self._clear_cache(ms1_tol, ms2_tol, score_thresh, min_match_peaks)

            self.speclib = None
            self._spec_hits = None
            self._bgc_hits = None
            self._rosetta_hits = []

        # next, try to load the cached rosetta_hits list. if parameters were invalidated above,
        # the file will have been deleted and this will fail
        logger.info("Trying to load cached Rosetta hits data")
        cached_rosetta_hits = load_pickled_data(self._nplinker, self._rhits_pickle_path)
        if cached_rosetta_hits is not None:
            logger.info(
                "Loaded cached Rosetta hits for dataset {} at {}".format(
                    self._dataset_id, self._rhits_pickle_path
                )
            )
            self._rosetta_hits = cached_rosetta_hits
            return self._rosetta_hits

        # if we get this far, it means regenerating some/all of the required data structures.
        #
        # create the _gnps2mibig and _mibig2gnps dicts if not already done
        if self._mibig2gnps is None or self._gnps2mibig is None:
            logger.info("Constructing GNPS/MiBIG dicts")
            self._load_csv(self._csv_path)

        self._init_bgc_hits(bgcs)

        # if we didn't find any BGC hits, no point in continuing
        if len(self._bgc_hits) == 0:
            logger.warning("Aborting Rosetta scoring data generation, no BGC hits were found!")
            # create an empty rosetta_hits.csv file
            self.export_to_csv(os.path.join(self._pickle_dir, "rosetta_hits.csv"))
            return self._rosetta_hits

        logger.info("No cached Rosetta hits data found")

        self._init_speclib(spectra, ms1_tol, ms2_tol, score_thresh, min_match_peaks)

        self._init_spec_hits(spectra, ms1_tol, ms2_tol, score_thresh, min_match_peaks)

        # finally construct the list of rosetta hits
        self._collect_rosetta_hits()

        # export cached data for future runs
        save_pickled_data(self._rosetta_hits, self._rhits_pickle_path)
        save_pickled_data(
            (Rosetta.PARAM_VERSION, ms1_tol, ms2_tol, score_thresh, min_match_peaks),
            self._params_pickle_path,
        )

        # automatically export CSV file containing hit data to <dataset>/rosetta
        # along with the pickled data
        self.export_to_csv(os.path.join(self._pickle_dir, "rosetta_hits.csv"))

        return self._rosetta_hits

    def _init_bgc_hits(self, bgcs):
        """Collect BGC hits. this is done first because the SpecLib generation below can take
        several minutes and is a waste of time if the knownclusterblast files required for
        the genomics data aren't available in the current dataset.
        """
        cached_bgc_hits = load_pickled_data(self._nplinker, self._bgchits_pickle_path)
        if cached_bgc_hits is not None and not self._ignore_genomic_cache:
            logger.info("Found pickled bgc_hits for dataset {}!".format(self._dataset_id))
            self._bgc_hits, self._mibig2bgc = cached_bgc_hits
        else:
            logger.info("Generating BGC hits")
            self._generate_bgc_hits(bgcs)

    def _init_spec_hits(self, spectra, ms1_tol, ms2_tol, score_thresh, min_match_peaks):
        spec_hits = load_pickled_data(self._nplinker, self._spechits_pickle_path)
        if spec_hits is not None:
            logger.info(
                "Found pickled spectral hits for dataset {} at {}".format(
                    self._dataset_id, self._spechits_pickle_path
                )
            )
            self._spec_hits = spec_hits

        if self._spec_hits is None:
            # no cached spectral hits available, generate (and cache)
            logger.info("Generating spectral hits")
            self._spec_hits = self._generate_spec_hits(
                spectra, ms1_tol, ms2_tol, score_thresh, min_match_peaks
            )

        logger.info(
            "SpecLib has {} spectra, {} hits".format(
                self.speclib.get_n_spec(), len(self._spec_hits)
            )
        )

    def _init_speclib(self, spectra, ms1_tol, ms2_tol, score_thresh, min_match_peaks):
        """Next is the metabolomic part. check if we have a pickled SpecLib object..."""
        speclib = load_pickled_data(self._nplinker, self._speclib_pickle_path)
        if speclib is not None:
            logger.info(
                "Found pickled SpecLib for dataset {} at {}!".format(
                    self._dataset_id, self._speclib_pickle_path
                )
            )
            self.speclib = speclib

        if self.speclib is None:
            # no cached speclib available, generate (and cache)
            logger.info("Generating SpecLib")
            self._generate_speclib(spectra)

    def _load_cached_params(self, ms1_tol, ms2_tol, score_thresh, min_match_peaks):
        """Check if cached parameters exist, and if so check they match the
        supplied ones. if not, need to regenerate any pickled data files.
        """
        params = load_pickled_data(self._nplinker, self._params_pickle_path)
        params_ok = False

        if params is not None:
            try:
                if params[0] != Rosetta.PARAM_VERSION:
                    logger.warning(
                        "Rosetta: pickled data version mismatch (old {}, new {})".format(
                            params[0], Rosetta.PARAM_VERSION
                        )
                    )
                else:
                    _version, _ms1_tol, _ms2_tol, _score_thresh, _min_match_peaks = params

                    if (
                        ms1_tol == _ms1_tol
                        and ms2_tol == _ms2_tol
                        and score_thresh == _score_thresh
                        and min_match_peaks == _min_match_peaks
                    ):
                        # params only valid if all of these match up
                        params_ok = True

            except Exception as e:
                logger.warning(f"Failed to parse pickled Rosetta parameters: {e}")

        return params_ok

    def _clear_cache(self, ms1_tol, ms2_tol, score_thresh, min_match_peaks):
        logger.info(
            "SpecLib parameters have been changed or do not exist, regenerating cached data files!"
        )
        logger.debug(
            "ms1_tol={:.3f}, ms2_tol={:.3f}, score_thresh={:.3f}, min_match_peaks={:d}".format(
                ms1_tol, ms2_tol, score_thresh, min_match_peaks
            )
        )
        for path in [
            self._bgchits_pickle_path,
            self._spechits_pickle_path,
            self._rhits_pickle_path,
            self._params_pickle_path,
            self._speclib_pickle_path,
            os.path.join(self._pickle_dir, "rosetta_hits.csv"),
        ]:
            if os.path.exists(path):
                os.unlink(path)

    def export_to_csv(self, filename):
        # convenience method for exporting a full set of rosetta hits to a CSV file
        with open(filename, "w", newline="") as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=",")

            csvwriter.writerow(
                [
                    "nplinker spectrum ID",
                    "spectrum ID",
                    "GNPS ID",
                    "spectral score",
                    "BGC ID",
                    "MiBIG BGC ID",
                    "BGC score",
                ]
            )
            for hit in self._rosetta_hits:
                csvwriter.writerow(
                    [
                        hit.spec.id,
                        hit.gnps_id,
                        hit.spec_match_score,
                        hit.bgc.id,
                        hit.mibig_id,
                        hit.bgc_match_score,
                    ]
                )
