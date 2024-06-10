import itertools
import logging
from nplinker.genomics.bgc import BGC
from nplinker.genomics.gcf import GCF
from nplinker.metabolomics import MolecularFamily
from nplinker.scoring.abc import ScoringBase
from nplinker.scoring.rosetta.rosetta import Rosetta
from .link_graph import LinkGraph
from .score import Score
from .scoring_method import ScoringMethod


logger = logging.getLogger(__name__)


class RosettaScoring(ScoringBase):
    name = ScoringMethod.ROSETTA.value
    ROSETTA_OBJ = None

    def __init__(self, npl):
        super().__init__(npl)
        self.bgc_to_gcf = True

        self.spec_score_cutoff = 0.0
        self.bgc_score_cutoff = 0.0

    @classmethod
    def setup(cls, npl):
        """Setup the Rosetta object and run the scoring algorithm.

        This method is only called once to setup the Rosetta object.
        """
        logger.info("RosettaScoring setup")
        cls.ROSETTA_OBJ = Rosetta(npl, ignore_genomic_cache=False)
        ms1_tol = Rosetta.DEF_MS1_TOL
        ms2_tol = Rosetta.DEF_MS2_TOL
        score_thresh = Rosetta.DEF_SCORE_THRESH
        min_match_peaks = Rosetta.DEF_MIN_MATCH_PEAKS

        ms1_tol, ms2_tol, score_thresh, min_match_peaks = RosettaScoring._init_from_config(
            npl.config
        )

        cls.ROSETTA_OBJ.run(npl.spectra, npl.bgcs, ms1_tol, ms2_tol, score_thresh, min_match_peaks)
        logger.info("RosettaScoring setup completed")

    @staticmethod
    def _init_from_config(config):
        """Allow overridding params via config file."""
        if "scoring" in config and "rosetta" in config["scoring"]:
            rc = config["scoring"]["rosetta"]
            ms1_tol = rc.get("ms1_tol", Rosetta.DEF_MS1_TOL)
            ms2_tol = rc.get("ms2_tol", Rosetta.DEF_MS2_TOL)
            score_thresh = rc.get("score_thresh", Rosetta.DEF_SCORE_THRESH)
            min_match_peaks = rc.get("min_match_peaks", Rosetta.DEF_MIN_MATCH_PEAKS)

        return ms1_tol, ms2_tol, score_thresh, min_match_peaks

    def _include_hit(self, hit):
        if (
            hit.spec_match_score < self.spec_score_cutoff
            or hit.bgc_match_score < self.bgc_score_cutoff
        ):
            return False

        return True

    def _insert_result_gen(self, results, src, hit):
        if src not in results:
            results[src] = {}
        # Rosetta can produce multiple "hits" per link
        if hit.spec in results[src]:
            original_data = results[src][hit.spec].data(self)
            results[src][hit.spec].set_data(self, original_data + [hit])
        else:
            results[src][hit.spec] = Score(name=self.name, value=[hit], parameter=self._params)

        return results

    def _insert_result_met(self, results, spec, target, hit):
        if spec not in results:
            results[spec] = {}
        # Rosetta can produce multiple "hits" per link
        if target in results[spec]:
            original_data = results[spec][target].data(self)
            results[spec][target].set_data(self, original_data + [hit])
        else:
            results[spec][target] = Score(name=self.name, value=[hit], parameter=self._params)

        return results

    def get_links(self, *objects, **parameters):
        # TODO: replace some attributes with parameters
        self._params = parameters

        self._validate_inputs(objects)

        if isinstance(objects[0], GCF):
            # assume user wants to use all BGCs from these GCFs
            bgcs = list(set(itertools.chain.from_iterable(x.bgcs for x in objects)))
            logger.info(
                "RosettaScoring got {} GCFs input, converted to {} BGCs".format(
                    len(objects), len(bgcs)
                )
            )
            objects = bgcs

        # list of RosettaHit objects which satisfy the current cutoffs
        ro_hits = list(
            filter(lambda hit: self._include_hit(hit), RosettaScoring.ROSETTA_OBJ._rosetta_hits)
        )

        # TODO this might need to be faster
        results = {}
        if isinstance(objects[0], BGC):
            results = self._collect_results_bgc(objects, ro_hits, results)
        else:  # Spectrum
            results = self._collect_results_spectra(objects, ro_hits, results)
        logger.debug(f"RosettaScoring found {len(results)} results")

        lg = LinkGraph()
        for src, links in results.items():
            for target, score in links.items():
                lg.add_link(src, target, score)
        return lg

    def _collect_results_spectra(self, objects, ro_hits, results):
        for spec in objects:
            for hit in ro_hits:
                if spec.id == hit.spec.id:
                    if not self.bgc_to_gcf:
                        # can use the BGC directly
                        results = self._insert_result_met(results, spec, hit.bgc, hit)
                    else:
                        # if we want the results to contain GCFs instead, need
                        # to iterate over the set of bgc.parents (most will only
                        # have one, hybrid BGCs will have more)
                        for gcf in hit.bgc.parents:
                            results = self._insert_result_met(results, spec, gcf, hit)
        return results

    def _collect_results_bgc(self, objects, ro_hits, results):
        for bgc in objects:
            for hit in ro_hits:
                if bgc.id == hit.bgc.id:
                    if not self.bgc_to_gcf:
                        # can use the BGC directly
                        results = self._insert_result_gen(results, bgc, hit)
                    else:
                        # if we want the results to contain GCFs instead, need
                        # to iterate over the set of bgc.parents (most will only
                        # have one, hybrid BGCs will have more)
                        for gcf in bgc.parents:
                            results = self._insert_result_gen(results, gcf, hit)
        return results

    def _validate_inputs(self, objects):
        """Enforce constraint that the list must contain a set of identically typed objects."""
        if not all(isinstance(x, type(objects[0])) for x in objects):
            raise Exception("RosettaScoring: uniformly-typed list of objects is required")

        if isinstance(objects[0], MolecularFamily):
            raise Exception("RosettaScoring requires input type Spectrum (found MolecularFamily)")

    def format_data(self, data):
        # TODO
        return f"{len(data)} hits"

    def sort(self, objects, reverse=True):
        # TODO
        return objects
