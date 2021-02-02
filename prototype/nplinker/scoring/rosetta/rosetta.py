import os
import csv

from ...logconfig import LogConfig
from ...parsers import mgf
from ...parsers.kcb import KCBParser

from ...pickler import load_pickled_data, save_pickled_data

from .rosetta_functions import fast_cosine
from .spec_lib import SpecLib

from ...genomics import MiBIGBGC

logger = LogConfig.getLogger(__file__)

class RosettaHit(object):

    def __init__(self, spec, gnps_id, mibig_id, bgc, spec_match_score, bgc_match_score):
        self.spec = spec
        self.gnps_id = gnps_id
        self.mibig_id = mibig_id
        self.bgc = bgc
        self.spec_match_score = spec_match_score
        self.bgc_match_score = bgc_match_score

    def __str__(self):
        return 'RosettaHit: {}<-->{} via ({} ({:.3f}), {} ({:.3f}))'.format(self.spec.spectrum_id, 
                                                                            self.bgc.name, 
                                                                            self.gnps_id, 
                                                                            self.spec_match_score,
                                                                            self.mibig_id,
                                                                            self.bgc_match_score)

    def __repr__(self):
        return str(self)

class Rosetta(object):

    DEF_MS1_TOL         = 100
    DEF_MS2_TOL         = 0.2
    DEF_SCORE_THRESH    = 0.5
    DEF_MIN_MATCH_PEAKS = 1

    PARAM_VERSION       = 1

    def __init__(self, nplinker, ignore_genomic_cache=False):
        self._nplinker = nplinker
        self._mgf_data = {}
        self._csv_data = {}
        self._mgf_path = os.path.join(nplinker.data_dir, 'matched_mibig_gnps_update.mgf')
        self._csv_path = os.path.join(nplinker.data_dir, 'matched_mibig_gnps_update.csv')
        self._data_path = nplinker.data_dir
        self._root_path = nplinker.root_dir
        self._dataset_id = nplinker.dataset_id
        self._ignore_genomic_cache = ignore_genomic_cache
        self._pickle_dir = os.path.join(nplinker.root_dir, 'rosetta')
        if not os.path.exists(self._pickle_dir):
            os.makedirs(self._pickle_dir, exist_ok=True)
        self._speclib_pickle_path = os.path.join(self._pickle_dir, 'SpecLib.pckl')
        self._bgchits_pickle_path = os.path.join(self._pickle_dir, 'bgc_hits.pckl')
        self._rhits_pickle_path = os.path.join(self._pickle_dir, 'RosettaHits.pckl')
        self._params_pickle_path = os.path.join(self._pickle_dir, 'RosettaParams.pckl')

        if not os.path.exists(self._mgf_path):
            logger.warning('Failed to load Rosetta data ({}), matching disabled'.format(self._mgf_path))
            return
        if not os.path.exists(self._csv_path):
            logger.warning('Failed to load Rosetta data ({}), matching disabled'.format(self._csv_path))
            return

        self._gnps2mibig = None
        self._mibig2gnps = None
        self._mibig2bgc = {}

        self.speclib = None
        self._rosetta_hits = []
        self._spec_hits = {}
        self._bgc_hits = {}

    @property
    def bgc_hits(self):
        return self._bgc_hits

    @property
    def spec_hits(self):
        return self._spec_hits

    def _load_csv(self, csv_path):
        logger.info('constructing rosetta dicts')

        self._gnps2mibig = {}
        self._mibig2gnps = {}

        with open(csv_path, 'r') as f:
            rdr = csv.reader(f, delimiter=',')
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
            hits = self.speclib.spectral_match(sp, ms2_tol=ms2_tol, min_match_peaks=min_match_peaks, ms1_tol=ms1_tol, score_thresh=score_thresh)
            if len(hits) > 0:
                spec_hits[sp] = hits
            if i % 100 == 0:
                logger.info('Searching for spectral hits {}/{}'.format(i, len(spectra)))

        return spec_hits

    def _load_speclib(self, spectra, ms1_tol, ms2_tol, score_thresh, min_match_peaks):
        logger.warning('No pickle SpecLib found, generating (this will take some time!)...')
        self.speclib = SpecLib(self._mgf_path)
        self.speclib._load_mgf()
        self.speclib.filter()

        logger.info('Finished generating SpecLib')
    
        # save_pickled_data((self.speclib, ms1_tol, ms2_tol, score_thresh, min_match_peaks), self._speclib_pickle_path)
        save_pickled_data(self.speclib, self._speclib_pickle_path)

    def _generate_bgc_hits(self, bgcs):
        self._bgc_hits = {}
        kcb_found = 0
        mibigs = 0
        errors = 0

        for i, bgc in enumerate(bgcs):
            # MiBIGBGC objects aren't created from local GenBank files, so don't
            # bother trying to lookup knownclusterblast results for them
            if isinstance(bgc, MiBIGBGC):
                mibigs += 1
                continue

            kcb_name = KCBParser.get_kcb_filename_from_bgc(bgc)
            if kcb_name is not None:
                try:
                    parser = KCBParser(kcb_name)
                    if len(parser.hits) > 0:
                        self._bgc_hits[bgc] = parser.hits
                    kcb_found += 1
                except Exception as e:
                    logger.warning(e)
                    errors += 1

            if i % 100 == 0:
                logger.info('Searching for BGC hits {}/{}'.format(i, len(bgcs)))

        # make reverse dict
        self._mibig2bgc = {}
        for bgc, hits in self._bgc_hits.items():
            for mibig_bgc_id in hits:
                if mibig_bgc_id not in self._mibig2bgc:
                    self._mibig2bgc[mibig_bgc_id] = set()
                self._mibig2bgc[mibig_bgc_id].add(bgc)

        logger.info('Completed, {} BGC hits found'.format(len(self._bgc_hits)))
        logger.info('GenBank files matched with knownclusterblast results: {}/{}'.format(kcb_found, len(bgcs) - mibigs))
        if errors > 0:
            logger.warning('Some knownclusterblast files could not be loaded, results may be incomplete')
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
                            self._rosetta_hits.append(RosettaHit(spec, gnps_id, mibig_id, bgc, score, bgc_score))
        logger.info('Found {} rosetta hits!'.format(len(self._rosetta_hits)))

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
                n_source_genes = len(hit[mibig_id]['all_bgc_genes'])
                n_mibig_genes = len(hit[mibig_id]['all_mibig_genes'])
                total_hit_identity = 0
                for hit_gene in hit[mibig_id]['individual_hits']:
                    identity_percent = hit_gene['identity_percent']
                    total_hit_identity += identity_percent / 100.0
                score = total_hit_identity / n_mibig_genes
                scores[mibig_id] = score
            processed[bgc] = scores
        return processed

    def _gather_kcb_filenames(self, antismash_path):
        pass

    
    def run(self, spectra, bgcs, ms1_tol, ms2_tol, score_thresh, min_match_peaks):
        # check if cached parameters exist, and if so check they match the 
        # supplied ones. if not, need to regenerate any pickled data files
        params = load_pickled_data(self._nplinker, self._params_pickle_path)
        params_ok = False

        if params is not None:
            try:
                if params[0] != Rosetta.PARAM_VERSION:
                    logger.warning('Rosetta: pickled data version mismatch (old {}, new {})'.format(params[0], Rosetta.PARAM_VERSION))
                else:
                    _version, _ms1_tol, _ms2_tol, _score_thresh, _min_match_peaks = params

                    if ms1_tol == _ms1_tol and ms2_tol == _ms2_tol and score_thresh == _score_thresh and min_match_peaks == _min_match_peaks:
                        # params only valid if all of these match up
                        params_ok = True

            except Exception as e:
                logger.warning('Failed to parse pickled Rosetta parameters: {}'.format(e))

        # if any parameters have been changed or version mismatch found, delete all cached files
        if not params_ok:
            logger.info('SpecLib parameters have been changed, regenerating cached data files!')
            logger.debug('ms1_tol={:.3f}, ms2_tol={:.3f}, score_thresh={:.3f}, min_match_peaks={:d}'.format(ms1_tol, ms2_tol, score_thresh, min_match_peaks))
            for path in [self._bgchits_pickle_path, self._rhits_pickle_path, self._params_pickle_path, self._speclib_pickle_path, os.path.join(self._pickle_dir, 'rosetta_hits.csv')]:
                if os.path.exists(path):
                    os.unlink(path)
            self.speclib = None
            self._rosetta_hits = []
            
        # next, try to load the cached rosetta_hits list. if parameters were invalidated above,
        # the file will have been deleted and this will fail
        logger.info('Trying to load cached Rosetta hits data')
        cached_rosetta_hits = load_pickled_data(self._nplinker, self._rhits_pickle_path)
        if cached_rosetta_hits is not None:
            logger.info('Loaded cached Rosetta hits for dataset {} at {}'.format(self._dataset_id, self._rhits_pickle_path))
            self._rosetta_hits = cached_rosetta_hits
            return self._rosetta_hits

        # if we get this far, it means regenerating some/all of the required data structures.
        # 
        # create the _gnps2mibig and _mibig2gnps dicts if not already done
        if self._mibig2gnps is None or self._gnps2mibig is None:
            logger.info('Constructing GNPS/MiBIG dicts')
            self._load_csv(self._csv_path)

        # collect BGC hits. this is done first because the SpecLib generation below can take 
        # several minutes and is a waste of time if the knownclusterblast files required for
        # the genomics data aren't available in the current dataset
        cached_bgc_hits = load_pickled_data(self._nplinker, self._bgchits_pickle_path)
        if cached_bgc_hits is not None and not self._ignore_genomic_cache:
            logger.info('Found pickled bgc_hits for dataset {}!'.format(self._dataset_id))
            self._bgc_hits, self._mibig2bgc = cached_bgc_hits
        else:
            logger.info('Generating BGC hits')
            self._generate_bgc_hits(bgcs)


        # if we didn't find any BGC hits, no point in continuing
        if len(self._bgc_hits) == 0:
            logger.warning('Aborting Rosetta scoring data generation, no BGC hits were found!')
            return self._rosetta_hits

        # next is the metabolomic part. check if we have a pickled SpecLib object...
        logger.info('No cached Rosetta hits data found')
        speclib = load_pickled_data(self._nplinker, self._speclib_pickle_path)
        if speclib is not None:
            logger.info('Found pickled SpecLib for dataset {} at {}!'.format(self._dataset_id, self._speclib_pickle_path))
            self.speclib = speclib
        
        if self.speclib is None:
            # no cached speclib available
            self._load_speclib(spectra, ms1_tol, ms2_tol, score_thresh, min_match_peaks)

        logger.info('Generating spectral hits')
        self._spec_hits = self._generate_spec_hits(spectra, ms1_tol, ms2_tol, score_thresh, min_match_peaks)

        logger.info('SpecLib has {} spectra, {} hits'.format(self.speclib.get_n_spec(), len(self._spec_hits)))

        # finally construct the list of rosetta hits
        self._collect_rosetta_hits()

        # export cached data for future runs
        save_pickled_data(self._rosetta_hits, self._rhits_pickle_path)
        save_pickled_data((Rosetta.PARAM_VERSION, ms1_tol, ms2_tol, score_thresh, min_match_peaks), self._params_pickle_path)

        # automatically export CSV file containing hit data to <dataset>/rosetta
        # along with the pickled data
        self.export_to_csv(os.path.join(self._pickle_dir, 'rosetta_hits.csv'))

        return self._rosetta_hits

    def export_to_csv(self, filename):
        # convenience method for exporting a full set of rosetta hits to a CSV file
        with open(filename, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=',')

            csvwriter.writerow(['nplinker spectrum ID', 'spectrum ID', 'GNPS ID', 'spectral score', 'nplinker BGC ID', 'BGC ID', 'MiBIG BGC ID', 'BGC score'])
            for hit in self._rosetta_hits:
                csvwriter.writerow([hit.spec.id, hit.spec.spectrum_id, hit.gnps_id, hit.spec_match_score, hit.bgc.id, hit.bgc.name, hit.mibig_id, hit.bgc_match_score])
