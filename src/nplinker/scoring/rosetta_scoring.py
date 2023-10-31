import itertools
from nplinker.genomics.bgc import BGC
from nplinker.genomics.gcf import GCF
from nplinker.logconfig import LogConfig
from nplinker.metabolomics import MolecularFamily
from nplinker.scoring.methods import ScoringMethod
from nplinker.scoring.object_link import ObjectLink
from nplinker.scoring.rosetta.rosetta import Rosetta


logger = LogConfig.getLogger(__name__)


class RosettaScoring(ScoringMethod):

    NAME = 'rosetta'
    ROSETTA_OBJ = None

    def __init__(self, npl):
        super().__init__(npl)
        self.bgc_to_gcf = True

        self.spec_score_cutoff = 0.0
        self.bgc_score_cutoff = 0.0

    @staticmethod
    def setup(npl):
        logger.info('RosettaScoring setup')
        RosettaScoring.ROSETTA_OBJ = Rosetta(npl, ignore_genomic_cache=False)
        ms1_tol = Rosetta.DEF_MS1_TOL
        ms2_tol = Rosetta.DEF_MS2_TOL
        score_thresh = Rosetta.DEF_SCORE_THRESH
        min_match_peaks = Rosetta.DEF_MIN_MATCH_PEAKS

        ms1_tol, ms2_tol, score_thresh, min_match_peaks = RosettaScoring._init_from_config(npl.config)

        RosettaScoring.ROSETTA_OBJ.run(npl.spectra, npl.bgcs, ms1_tol, ms2_tol,
                                       score_thresh, min_match_peaks)
        logger.info('RosettaScoring setup completed')

    @staticmethod
    def _init_from_config(config):
        """ allow overridding params via config file
        """
        if 'scoring' in config and 'rosetta' in config['scoring']:
            rc = config['scoring']['rosetta']
            ms1_tol = rc.get('ms1_tol', Rosetta.DEF_MS1_TOL)
            ms2_tol = rc.get('ms2_tol', Rosetta.DEF_MS2_TOL)
            score_thresh = rc.get('score_thresh', Rosetta.DEF_SCORE_THRESH)
            min_match_peaks = rc.get('min_match_peaks',
                                     Rosetta.DEF_MIN_MATCH_PEAKS)

        return ms1_tol,ms2_tol,score_thresh,min_match_peaks

    def _include_hit(self, hit):
        if hit.spec_match_score < self.spec_score_cutoff or hit.bgc_match_score < self.bgc_score_cutoff:
            return False

        return True

    def _insert_result_gen(self, results, src, hit):
        if src not in results:
            results[src] = {}
        # Rosetta can produce multiple "hits" per link, need to
        # ensure the ObjectLink contains all the RosettaHit objects
        # in these cases
        if hit.spec in results[src]:
            original_data = results[src][hit.spec].data(self)
            results[src][hit.spec].set_data(self, original_data + [hit])
        else:
            results[src][hit.spec] = ObjectLink(src,
                                                hit.spec,
                                                self,
                                                data=[hit])

        return results

    def _insert_result_met(self, results, spec, target, hit):
        if spec not in results:
            results[spec] = {}
        # Rosetta can produce multiple "hits" per link, need to
        # ensure the ObjectLink contains all the RosettaHit objects
        # in these cases
        if target in results[spec]:
            original_data = results[spec][target].data(self)
            results[spec][target].set_data(self, original_data + [hit])
        else:
            results[spec][target] = ObjectLink(spec, target, self, data=[hit])

        return results

    def get_links(self, objects, link_collection):
        self._validate_inputs(objects)

        if isinstance(objects[0], GCF):
            # assume user wants to use all BGCs from these GCFs
            bgcs = list(
                set(itertools.chain.from_iterable(x.bgcs for x in objects)))
            logger.info(
                'RosettaScoring got {} GCFs input, converted to {} BGCs'.
                format(len(objects), len(bgcs)))
            objects = bgcs

        # list of RosettaHit objects which satisfy the current cutoffs
        ro_hits = list(
            filter(lambda hit: self._include_hit(hit),
                   RosettaScoring.ROSETTA_OBJ._rosetta_hits))

        # TODO this might need to be faster
        results = {}
        if isinstance(objects[0], BGC):
            results = self._collect_results_bgc(objects, ro_hits, results)
        else:  # Spectrum
            results = self._collect_results_spectra(objects, ro_hits, results)

        link_collection._add_links_from_method(self, results)
        logger.debug(f'RosettaScoring found {len(results)} results')
        return link_collection

    def _collect_results_spectra(self, objects, ro_hits, results):
        for spec in objects:
            for hit in ro_hits:
                if spec.id == hit.spec.id:
                    if not self.bgc_to_gcf:
                            # can use the BGC directly
                        results = self._insert_result_met(
                                results, spec, hit.bgc, hit)
                    else:
                            # if we want the results to contain GCFs instead, need
                            # to iterate over the set of bgc.parents (most will only
                            # have one, hybrid BGCs will have more)
                        for gcf in hit.bgc.parents:
                            results = self._insert_result_met(
                                    results, spec, gcf, hit)
        return results

    def _collect_results_bgc(self, objects, ro_hits, results):
        for bgc in objects:
            for hit in ro_hits:
                if bgc.bgc_id == hit.bgc.bgc_id:
                    if not self.bgc_to_gcf:
                            # can use the BGC directly
                        results = self._insert_result_gen(
                                results, bgc, hit)
                    else:
                            # if we want the results to contain GCFs instead, need
                            # to iterate over the set of bgc.parents (most will only
                            # have one, hybrid BGCs will have more)
                        for gcf in bgc.parents:
                            results = self._insert_result_gen(
                                    results, gcf, hit)
        return results

    def _validate_inputs(self, objects):
        """ enforce constraint that the list must contain a set of identically typed objects
        """
        if not all(isinstance(x, type(objects[0])) for x in objects):
            raise Exception(
                'RosettaScoring: uniformly-typed list of objects is required')

        if isinstance(objects[0], MolecularFamily):
            raise Exception(
                'RosettaScoring requires input type Spectrum (found MolecularFamily)'
            )

    def format_data(self, data):
        # TODO
        return f'{len(data)} hits'

    def sort(self, objects, reverse=True):
        # TODO
        return objects
