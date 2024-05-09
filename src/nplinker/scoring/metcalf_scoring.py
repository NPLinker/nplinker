from __future__ import annotations
import logging
import os
from typing import TYPE_CHECKING
import numpy as np
import pandas as pd
from nplinker.defaults import OUTPUT_DEFAULT_PATH
from nplinker.genomics import GCF
from nplinker.metabolomics import MolecularFamily
from nplinker.metabolomics import Spectrum
from nplinker.pickler import load_pickled_data
from nplinker.pickler import save_pickled_data
from .linking import LINK_TYPES
from .linking import DataLinks
from .linking import LinkFinder
from .linking import isinstance_all
from .methods import ScoringMethod
from .object_link import ObjectLink


if TYPE_CHECKING:
    from ..nplinker import NPLinker
    from . import LinkCollection

logger = logging.getLogger(__name__)


class MetcalfScoring(ScoringMethod):
    """Metcalf scoring method.

    Attributes:
        DATALINKS: The DataLinks object to use for scoring.
        LINKFINDER: The LinkFinder object to use for scoring.
        NAME: The name of the scoring method. This is set to 'metcalf'.
        CACHE: The name of the cache file to use for storing the MetcalfScoring.
    """

    DATALINKS = None
    LINKFINDER = None
    NAME = "metcalf"
    CACHE = "cache_metcalf_scoring.pckl"

    def __init__(self, npl: NPLinker) -> None:
        """Create a MetcalfScoring object.

        Args:
            npl: The NPLinker object to use for scoring.

        Attributes:
            cutoff: The cutoff value to use for scoring. Scores below
                this value will be discarded. Defaults to 1.0.
            standardised: Whether to use standardised scores. Defaults
                to True.
            name: The name of the scoring method. It's set to a fixed value
                'metcalf'.
        """
        super().__init__(npl)
        self.cutoff = 1.0
        self.standardised = True

    # TODO CG: not sure why using staticmethod here. Check later and refactor if possible
    # TODO CG: refactor this method and extract code for cache file to a separate method
    @staticmethod
    def setup(npl: NPLinker):
        """Setup the MetcalfScoring object.

        DataLinks and LinkFinder objects are created and cached for later use.
        """
        logger.info(
            "MetcalfScoring.setup (bgcs={}, gcfs={}, spectra={}, molfams={}, strains={})".format(
                len(npl.bgcs), len(npl.gcfs), len(npl.spectra), len(npl.molfams), len(npl.strains)
            )
        )

        OUTPUT_DEFAULT_PATH.mkdir(exist_ok=True)
        cache_file = OUTPUT_DEFAULT_PATH / MetcalfScoring.CACHE

        # the metcalf preprocessing can take a long time for large datasets, so it's
        # better to cache as the data won't change unless the number of objects does
        dataset_counts = [
            len(npl.bgcs),
            len(npl.gcfs),
            len(npl.spectra),
            len(npl.molfams),
            len(npl.strains),
        ]
        datalinks, linkfinder = None, None
        if os.path.exists(cache_file):
            logger.info("MetcalfScoring.setup loading cached data")
            cache_data = load_pickled_data(npl, cache_file)
            cache_ok = True
            if cache_data is not None:
                (counts, datalinks, linkfinder) = cache_data
                # need to invalidate this if dataset appears to have changed
                for i in range(len(counts)):
                    if counts[i] != dataset_counts[i]:
                        logger.info("MetcalfScoring.setup invalidating cached data!")
                        cache_ok = False
                        break

            if cache_ok:
                MetcalfScoring.DATALINKS = datalinks
                MetcalfScoring.LINKFINDER = linkfinder

        if MetcalfScoring.DATALINKS is None:
            logger.info("MetcalfScoring.setup preprocessing dataset (this may take some time)")
            MetcalfScoring.DATALINKS = DataLinks(npl.gcfs, npl.spectra, npl.molfams, npl.strains)
            MetcalfScoring.LINKFINDER = LinkFinder()
            MetcalfScoring.LINKFINDER.calc_score(MetcalfScoring.DATALINKS, link_type=LINK_TYPES[0])
            MetcalfScoring.LINKFINDER.calc_score(MetcalfScoring.DATALINKS, link_type=LINK_TYPES[1])
            logger.info("MetcalfScoring.setup caching results")
            save_pickled_data(
                (dataset_counts, MetcalfScoring.DATALINKS, MetcalfScoring.LINKFINDER), cache_file
            )

        logger.info("MetcalfScoring.setup completed")

    # TODO CG: is it needed? remove it if not
    @property
    def datalinks(self) -> DataLinks | None:
        """Get the DataLinks object used for scoring."""
        return MetcalfScoring.DATALINKS

    def get_links(
        self, *objects: GCF | Spectrum | MolecularFamily, link_collection: LinkCollection
    ) -> LinkCollection:
        """Get links for the given objects and add them to the given LinkCollection.

        The given objects are treated as input or source objects, which must
        be GCF, Spectrum or MolecularFamily objects.

        Args:
            objects: The objects to get links for. Must be GCF, Spectrum
                or MolecularFamily objects.
            link_collection: The LinkCollection object to add the links to.

        Returns:
            The LinkCollection object with the new links added.

        Raises:
            ValueError: If the input objects are empty.
            TypeError: If the input objects are not of the correct type.
            ValueError: If LinkFinder instance has not been created
                (MetcalfScoring object has not been setup).
        """
        if len(objects) == 0:
            raise ValueError("Empty input objects.")

        if isinstance_all(*objects, objtype=GCF):
            obj_type = "gcf"
        elif isinstance_all(*objects, objtype=Spectrum):
            obj_type = "spec"
        elif isinstance_all(*objects, objtype=MolecularFamily):
            obj_type = "mf"
        else:
            types = [type(i) for i in objects]
            raise TypeError(
                f"Invalid type {set(types)}. Input objects must be GCF, Spectrum or MolecularFamily objects."
            )

        if self.LINKFINDER is None:
            raise ValueError(
                ("LinkFinder object not found. Have you called `MetcalfScoring.setup(npl)`?")
            )

        logger.info(f"MetcalfScoring: standardised = {self.standardised}")
        if not self.standardised:
            scores_list = self.LINKFINDER.get_links(*objects, score_cutoff=self.cutoff)
        # TODO CG: verify the logics of standardised score and add unit tests
        else:
            # use negative infinity as the score cutoff to ensure we get all links
            # the self.cutoff will be applied later in the postprocessing step
            scores_list = self.LINKFINDER.get_links(*objects, score_cutoff=np.NINF)
            if obj_type == "gcf":
                scores_list = self._calc_standardised_score_gen(self.LINKFINDER, scores_list)
            else:
                scores_list = self._calc_standardised_score_met(self.LINKFINDER, scores_list)

        link_scores: dict[
            GCF | Spectrum | MolecularFamily, dict[GCF | Spectrum | MolecularFamily, ObjectLink]
        ] = {}
        if obj_type == "gcf":
            logger.info(
                f"MetcalfScoring: input_type=GCF, result_type=Spec/MolFam, "
                f"#inputs={len(objects)}."
            )
            for scores in scores_list:
                # when no links found
                if scores.shape[1] == 0:
                    logger.info(f'MetcalfScoring: found no "{scores.name}" links')
                else:
                    # when links found
                    for col_index in range(scores.shape[1]):
                        gcf = self.npl.lookup_gcf(scores.loc["source", col_index])
                        if scores.name == LINK_TYPES[0]:
                            met = self.npl.lookup_spectrum(scores.loc["target", col_index])
                        else:
                            met = self.npl.lookup_mf(scores.loc["target", col_index])
                        if gcf not in link_scores:
                            link_scores[gcf] = {}
                        # TODO CG: use id instead of object for gcf, met and self?
                        link_scores[gcf][met] = ObjectLink(
                            gcf, met, self, scores.loc["score", col_index]
                        )
                    logger.info(f"MetcalfScoring: found {len(link_scores)} {scores.name} links.")
        else:
            logger.info(
                f"MetcalfScoring: input_type=Spec/MolFam, result_type=GCF, "
                f"#inputs={len(objects)}."
            )
            scores = scores_list[0]
            # when no links found
            if scores.shape[1] == 0:
                logger.info(f'MetcalfScoring: found no links "{scores.name}" for input objects')
            else:
                for col_index in range(scores.shape[1]):
                    gcf = self.npl.lookup_gcf(scores.loc["target", col_index])
                    if scores.name == LINK_TYPES[0]:
                        met = self.npl.lookup_spectrum(scores.loc["source", col_index])
                    else:
                        met = self.npl.lookup_mf(scores.loc["source", col_index])
                    if met not in link_scores:
                        link_scores[met] = {}
                    link_scores[met][gcf] = ObjectLink(
                        met, gcf, self, scores.loc["score", col_index]
                    )
                logger.info(f"MetcalfScoring: found {len(link_scores)} {scores.name} links.")

        link_collection._add_links_from_method(self, link_scores)
        logger.info("MetcalfScoring: completed")
        return link_collection

    def _calc_standardised_score_met(
        self, linkfinder: LinkFinder, results: list
    ) -> list[pd.DataFrame]:
        if linkfinder.metcalf_mean is None or linkfinder.metcalf_std is None:
            raise ValueError(
                "Metcalf mean and std not found. Have you called `MetcalfScoring.setup(npl)`?"
            )
        logger.info("Calculating standardised Metcalf scores (met input)")
        raw_score = results[0]
        z_scores = []
        for col_index in range(raw_score.shape[1]):
            gcf = self.npl.lookup_gcf(raw_score.loc["target", col_index])
            if raw_score.name == LINK_TYPES[0]:
                met = self.npl.lookup_spectrum(raw_score.at["source", col_index])
            else:
                met = self.npl.lookup_mf(raw_score.at["source", col_index])

            num_gcf_strains = len(gcf.strains)
            num_met_strains = len(met.strains)
            mean = linkfinder.metcalf_mean[num_met_strains][num_gcf_strains]
            sqrt = linkfinder.metcalf_std[num_met_strains][num_gcf_strains]
            z_score = (raw_score.at["score", col_index] - mean) / sqrt
            z_scores.append(z_score)

        z_scores = np.array(z_scores)
        mask = z_scores >= self.cutoff

        scores_df = pd.DataFrame(
            [
                raw_score.loc["source"].values[mask],
                raw_score.loc["target"].values[mask],
                z_scores[mask],
            ],
            index=raw_score.index,
        )
        scores_df.name = raw_score.name

        return [scores_df]

    def _calc_standardised_score_gen(
        self, linkfinder: LinkFinder, results: list
    ) -> list[pd.DataFrame]:
        if linkfinder.metcalf_mean is None or linkfinder.metcalf_std is None:
            raise ValueError(
                "Metcalf mean and std not found. Have you called `MetcalfScoring.setup(npl)`?"
            )
        logger.info("Calculating standardised Metcalf scores (gen input)")
        postprocessed_scores = []
        for raw_score in results:
            z_scores = []
            for col_index in range(raw_score.shape[1]):
                gcf = self.npl.lookup_gcf(raw_score.loc["source", col_index])
                if raw_score.name == LINK_TYPES[0]:
                    met = self.npl.lookup_spectrum(raw_score.at["target", col_index])
                else:
                    met = self.npl.lookup_mf(raw_score.at["target", col_index])

                num_gcf_strains = len(gcf.strains)
                num_met_strains = len(met.strains)
                mean = linkfinder.metcalf_mean[num_met_strains][num_gcf_strains]
                sqrt = linkfinder.metcalf_std[num_met_strains][num_gcf_strains]
                z_score = (raw_score.at["score", col_index] - mean) / sqrt
                z_scores.append(z_score)

            z_scores = np.array(z_scores)
            mask = z_scores >= self.cutoff

            scores_df = pd.DataFrame(
                [
                    raw_score.loc["source"].values[mask],
                    raw_score.loc["target"].values[mask],
                    z_scores[mask],
                ],
                index=raw_score.index,
            )
            scores_df.name = raw_score.name
            postprocessed_scores.append(scores_df)

        return postprocessed_scores

    # TODO CG: refactor this method
    def format_data(self, data):
        """Format the data for display."""
        # for metcalf the data will just be a floating point value (i.e. the score)
        return f"{data:.4f}"

    # TODO CG: refactor this method
    def sort(self, objects, reverse=True):
        """Sort the objects based on the score."""
        # sort based on score
        return sorted(objects, key=lambda objlink: objlink[self], reverse=reverse)
