import os
import csv
import pickle

from ...logconfig import LogConfig
from ...parsers import mgf
from ...parsers.kcb import KCBParser

from .rosetta_functions import fast_cosine
from .spec_lib import SpecLib

logger = LogConfig.getLogger(__file__)

class RosettaHit(object):

    def __init__(self, spec, gnps_id, mibig_id, bgc):
        self.spec = spec
        self.gnps_id = gnps_id
        self.mibig_id = mibig_id
        self.bgc = bgc

    def __str__(self):
        return 'RosettaHit: {}<-->{} via ({}, {})'.format(self.spec.spectrum_id, self.bgc.name, self.gnps_id, self.mibig_id)

    def __repr__(self):
        return str(self)

class Rosetta(object):

    def __init__(self, data_path, root_path, dataset_id):
        self._mgf_data = {}
        self._csv_data = {}
        self._mgf_path = os.path.join(data_path, 'matched_mibig_gnps_update.mgf')
        self._csv_path = os.path.join(data_path, 'matched_mibig_gnps_update.csv')
        self._data_path = data_path
        self._root_path = root_path
        self._dataset_id = dataset_id
        self._pickle_dir = os.path.join(root_path, 'rosetta')
        if not os.path.exists(self._pickle_dir):
            os.makedirs(self._pickle_dir, exist_ok=True)
        self._speclib_pickle_path = os.path.join(self._pickle_dir, dataset_id + '_SpecLib.pckl')
        self._bgchits_pickle_path = os.path.join(self._pickle_dir, dataset_id + '_bgc_hits.pckl')

        if not os.path.exists(self._mgf_path):
            logger.warning('Failed to load Rosetta data ({}), matching disabled'.format(self._mgf_path))
            return
        if not os.path.exists(self._csv_path):
            logger.warning('Failed to load Rosetta data ({}), matching disabled'.format(self._csv_path))
            return

        self._gnps2mibig = None
        self._mibig2gnps = None
        self._mibig2bgc = {}

        self._rosetta_hits = []
        self._spec_hits = {}
        self._bgc_hits = {}

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

    def _generate_spec_hits(self, spectra, ms1_tol, score_thresh):
        spec_hits = {}
        for i, sp in enumerate(spectra):
            hits = self.speclib.spectral_match(sp, ms1_tol=ms1_tol, score_thresh=score_thresh)
            if len(hits) > 0:
                spec_hits[sp] = hits
            if i % 100 == 0:
                logger.info('Searching for spectral hits {}/{}'.format(i, len(spectra)))

        return spec_hits

    def _load_speclib(self, spectra):
        logger.warning('No pickle SpecLib found, generating (this will take some time!)...')
        self.speclib = SpecLib(self._mgf_path)
        self.speclib._load_mgf()
        self.speclib.filter()

        logger.info('Finished generating SpecLib')
    
        with open(self._speclib_pickle_path, 'wb') as f:
            pickle.dump(self.speclib, f)

    def _generate_bgc_hits(self, bgcs):
        self._bgc_hits = {}
        for i, bgc in enumerate(bgcs):
            bgc_file = bgc.antismash_file
            kcb_name = KCBParser.get_kcb_filename_from_bgc(bgc)
            if kcb_name is not None:
                parser = KCBParser(kcb_name)
                if len(parser.hits) > 0:
                    self._bgc_hits[bgc] = parser.hits

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
        with open(self._bgchits_pickle_path, 'wb') as f:
            pickle.dump((self._bgc_hits, self._mibig2bgc), f)

    def _collect_rosetta_hits(self):
        self._rosetta_hits = []
        for spec, data in self._spec_hits.items():
            for gnps_id, score in data:
                for mibig_id in self._gnps2mibig[gnps_id]:
                    if mibig_id in self._mibig2bgc:
                        for bgc in self._mibig2bgc[mibig_id]:
                            self._rosetta_hits.append(RosettaHit(spec, gnps_id, mibig_id, bgc))
        logger.info('Found {} rosetta hits!'.format(len(self._rosetta_hits)))

    def generate_bgc_summary_scores(self):
        # TODO
        # this method is currently broken due to missing 'all_mibig_genes'
        # in knownclusterblast parser output...

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
                n_source_genes = len(hit[mibig_id][0]['all_bgc_genes'])
                n_mibig_genes = len(hit[mibig_id][0]['all_mibig_genes'])
                total_hit_identity = 0
                for hit_gene in hit[mibig_id]:
                    identity_percent = hit_gene['identity_percent']
                    total_hit_identity += identity_percent / 100.0
                score = total_hit_identity / n_mibig_genes
                scores[mibig_id] = score
            processed[hit] = scores
        return scores

    def run(self, spectra, bgcs, ms1_tol=100, score_thresh=0.5):
        # check if we have a pickled SpecLib object first...
        if os.path.exists(self._speclib_pickle_path):
            logger.info('Found pickled SpecLib for dataset {}!'.format(self._dataset_id))
            with open(self._speclib_pickle_path, 'rb') as f:
                self.speclib = pickle.load(f)
        else:
            self._load_speclib(spectra)

        self._spec_hits = self._generate_spec_hits(spectra, ms1_tol, score_thresh)

        logger.info('SpecLib has {} spectra, {} hits'.format(self.speclib.get_n_spec(), len(self._spec_hits)))

        # create the _gnps2mibig and _mibig2gnps dicts if not already done
        if self._mibig2gnps is None or self._gnps2mibig is None:
            self._load_csv(self._csv_path)

        # collect bgc hits
        if os.path.exists(self._bgchits_pickle_path):
            logger.info('Found pickled bgc_hits for dataset {}!'.format(self._dataset_id))
            with open(self._bgchits_pickle_path, 'rb') as f:
                self._bgc_hits, self._mibig2bgc = pickle.load(f)
        else:
            self._generate_bgc_hits(bgcs)

        # finally construct the list of rosetta hits
        self._collect_rosetta_hits()

        return self._rosetta_hits
