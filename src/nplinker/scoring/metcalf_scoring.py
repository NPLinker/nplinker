from __future__ import annotations
import logging
import os
from typing import TYPE_CHECKING
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from nplinker.genomics import GCF
from nplinker.metabolomics import MolecularFamily
from nplinker.metabolomics import Spectrum
from nplinker.pickler import load_pickled_data
from nplinker.pickler import save_pickled_data
from .abc import ScoringBase
from .linking import LINK_TYPES
from .linking import DataLinks
from .linking import isinstance_all
from .object_link import ObjectLink


if TYPE_CHECKING:
    from ..nplinker import NPLinker
    from . import LinkCollection

logger = logging.getLogger(__name__)


class MetcalfScoring(ScoringBase):
    """Metcalf scoring method.

    Attributes:
        name: The name of this scoring method, set to a fixed value `metcalf`.
        DATALINKS: The DataLinks object to use for scoring.
        CACHE: The name of the cache file to use for storing the MetcalfScoring.

        raw_score_spec_gcf: The raw Metcalf scores for spectrum-GCF links.
        raw_score_mf_gcf: The raw Metcalf scores for molecular family-GCF links.
        metcalf_mean: The mean value used for standardising Metcalf scores.
        metcalf_std: The standard deviation value used for standardising Metcalf scores.
    """

    name = "metcalf"
    DATALINKS = None
    CACHE = "cache_metcalf_scoring.pckl"

    raw_score_spec_gcf = pd.DataFrame()
    raw_score_mf_gcf = pd.DataFrame()
    metcalf_mean = None
    metcalf_std = None

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

    # TODO CG: refactor this method and extract code for cache file to a separate method
    @classmethod
    def setup(cls, npl: NPLinker):
        """Setup the DataLinks object.

        This method is only called once to setup the DataLinks object.
        """
        logger.info(
            "MetcalfScoring.setup (bgcs={}, gcfs={}, spectra={}, molfams={}, strains={})".format(
                len(npl.bgcs), len(npl.gcfs), len(npl.spectra), len(npl.molfams), len(npl.strains)
            )
        )

        cache_file = npl.output_dir / cls.CACHE

        # the metcalf preprocessing can take a long time for large datasets, so it's
        # better to cache as the data won't change unless the number of objects does
        dataset_counts = [
            len(npl.bgcs),
            len(npl.gcfs),
            len(npl.spectra),
            len(npl.molfams),
            len(npl.strains),
        ]
        datalinks = None
        if os.path.exists(cache_file):
            logger.info("MetcalfScoring.setup loading cached data")
            cache_data = load_pickled_data(npl, cache_file)
            cache_ok = True
            if cache_data is not None:
                (counts, datalinks) = cache_data
                # need to invalidate this if dataset appears to have changed
                for i in range(len(counts)):
                    if counts[i] != dataset_counts[i]:
                        logger.info("MetcalfScoring.setup invalidating cached data!")
                        cache_ok = False
                        break

            if cache_ok:
                cls.DATALINKS = datalinks

        if cls.DATALINKS is None:
            logger.info("MetcalfScoring.setup preprocessing dataset (this may take some time)")
            cls.DATALINKS = DataLinks(npl.gcfs, npl.spectra, npl.molfams, npl.strains)
            cls.calc_score(cls.DATALINKS, link_type=LINK_TYPES[0])
            cls.calc_score(cls.DATALINKS, link_type=LINK_TYPES[1])
            logger.info("MetcalfScoring.setup caching results")
            # TODO: save the score values 2024-05-29
            save_pickled_data((dataset_counts, cls.DATALINKS), cache_file)

        logger.info("MetcalfScoring.setup completed")

    @classmethod
    def calc_score(
        cls,
        data_links: DataLinks,
        link_type: str = "spec-gcf",
        scoring_weights: tuple[int, int, int, int] = (10, -10, 0, 1),
    ) -> None:
        """Calculate Metcalf scores.

        This method calculates the `raw_score_spec_gcf`, `raw_score_mf_gcf`, `metcalf_mean`, and
        `metcalf_std` attributes.

        Args:
            data_links: The DataLinks object to use for scoring.
            link_type: The type of link to score. Must be 'spec-gcf' or
                'mf-gcf'. Defaults to 'spec-gcf'.
            scoring_weights: The weights to
                use for Metcalf scoring. The weights are applied to
                '(met_gcf, met_not_gcf, gcf_not_met, not_met_not_gcf)'.
                Defaults to (10, -10, 0, 1).

        Raises:
            ValueError: If an invalid link type is provided.
        """
        if link_type not in LINK_TYPES:
            raise ValueError(f"Invalid link type: {link_type}. Must be one of {LINK_TYPES}")

        if link_type == "spec-gcf":
            cls.raw_score_spec_gcf = (
                data_links.cooccurrence_spec_gcf * scoring_weights[0]
                + data_links.cooccurrence_spec_notgcf * scoring_weights[1]
                + data_links.cooccurrence_notspec_gcf * scoring_weights[2]
                + data_links.cooccurrence_notspec_notgcf * scoring_weights[3]
            )

        if link_type == "mf-gcf":
            cls.raw_score_mf_gcf = (
                data_links.cooccurrence_mf_gcf * scoring_weights[0]
                + data_links.cooccurrence_mf_notgcf * scoring_weights[1]
                + data_links.cooccurrence_notmf_gcf * scoring_weights[2]
                + data_links.cooccurrence_notmf_notgcf * scoring_weights[3]
            )

        if cls.metcalf_mean is None or cls.metcalf_std is None:
            n_strains = data_links.occurrence_gcf_strain.shape[1]
            cls.metcalf_mean, cls.metcalf_std = cls._calc_mean_std(n_strains, scoring_weights)

    # TODO CG: is it needed? remove it if not
    @property
    def datalinks(self) -> DataLinks | None:
        """Get the DataLinks object used for scoring."""
        return self.DATALINKS

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

        logger.info(f"MetcalfScoring: standardised = {self.standardised}")
        if not self.standardised:
            scores_list = self._get_links(*objects, score_cutoff=self.cutoff)
        # TODO CG: verify the logics of standardised score and add unit tests
        else:
            # use negative infinity as the score cutoff to ensure we get all links
            # the self.cutoff will be applied later in the postprocessing step
            scores_list = self._get_links(*objects, score_cutoff=np.NINF)
            if obj_type == "gcf":
                scores_list = self._calc_standardised_score_gen(scores_list)
            else:
                scores_list = self._calc_standardised_score_met(scores_list)

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

    @staticmethod
    def _calc_mean_std(
        n_strains: int, scoring_weights: tuple[int, int, int, int]
    ) -> tuple[np.ndarray, np.ndarray]:
        sz = (n_strains + 1, n_strains + 1)
        mean = np.zeros(sz)
        variance = np.zeros(sz)
        for n in range(n_strains + 1):
            for m in range(n_strains + 1):
                max_overlap = min(n, m)
                min_overlap = max(0, n + m - n_strains)
                expected_value = 0
                expected_sq = 0
                for o in range(min_overlap, max_overlap + 1):
                    o_prob = hypergeom.pmf(o, n_strains, n, m)
                    # compute metcalf for n strains in type 1 and m in gcf
                    score = o * scoring_weights[0]
                    score += scoring_weights[1] * (n - o)
                    score += scoring_weights[2] * (m - o)
                    score += scoring_weights[3] * (n_strains - (n + m - o))
                    expected_value += o_prob * score
                    expected_sq += o_prob * (score**2)
                mean[n, m] = expected_value
                expected_sq = expected_sq - expected_value**2
                if expected_sq < 1e-09:
                    expected_sq = 1
                variance[n, m] = expected_sq
        return mean, np.sqrt(variance)

    def _get_links(
        self,
        *objects: tuple[GCF, ...] | tuple[Spectrum, ...] | tuple[MolecularFamily, ...],
        score_cutoff: float = 0.5,
    ) -> list[pd.DataFrame]:
        """Get links and scores for given objects.

        Args:
            objects: A list of GCF, Spectrum or MolecularFamily objects
                and all objects must be of the same type.
            score_cutoff: Minimum score to consider a link (≥score_cutoff).
                Default is 0.5.

        Returns:
            List of data frames containing the ids of the linked objects
                and the score. The data frame has index names of
                'source', 'target' and 'score':

                - the 'source' row contains the ids of the input/source objects,
                - the 'target' row contains the ids of the target objects,
                - the 'score' row contains the scores.

        Raises:
            ValueError: If input objects are empty.
            TypeError: If input objects are not GCF, Spectrum or MolecularFamily objects.
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

        links = []
        if obj_type == "gcf":
            # TODO CG: the hint and mypy warnings will be gone after renaming all
            # string ids to `.id`
            obj_ids = [gcf.gcf_id for gcf in objects]
            # spec-gcf
            scores = self.raw_score_spec_gcf.loc[:, obj_ids]
            df = self._get_scores_source_gcf(scores, score_cutoff)
            df.name = LINK_TYPES[0]
            links.append(df)
            # mf-gcf
            scores = self.raw_score_mf_gcf.loc[:, obj_ids]
            df = self._get_scores_source_gcf(scores, score_cutoff)
            df.name = LINK_TYPES[1]
            links.append(df)

        if obj_type == "spec":
            obj_ids = [spec.spectrum_id for spec in objects]
            scores = self.raw_score_spec_gcf.loc[obj_ids, :]
            df = self._get_scores_source_met(scores, score_cutoff)
            df.name = LINK_TYPES[0]
            links.append(df)

        if obj_type == "mf":
            obj_ids = [mf.family_id for mf in objects]
            scores = self.raw_score_mf_gcf.loc[obj_ids, :]
            df = self._get_scores_source_met(scores, score_cutoff)
            df.name = LINK_TYPES[1]
            links.append(df)
        return links

    @staticmethod
    def _get_scores_source_gcf(scores: pd.DataFrame, score_cutoff: float) -> pd.DataFrame:
        row_indexes, col_indexes = np.where(scores >= score_cutoff)
        src_obj_ids = scores.columns[col_indexes].to_list()
        target_obj_ids = scores.index[row_indexes].to_list()
        scores_candidate = scores.values[row_indexes, col_indexes].tolist()
        return pd.DataFrame(
            [src_obj_ids, target_obj_ids, scores_candidate], index=["source", "target", "score"]
        )

    @staticmethod
    def _get_scores_source_met(scores: pd.DataFrame, score_cutoff: float) -> pd.DataFrame:
        row_indexes, col_indexes = np.where(scores >= score_cutoff)
        src_obj_ids = scores.index[row_indexes].to_list()
        target_obj_ids = scores.columns[col_indexes].to_list()
        scores_candidate = scores.values[row_indexes, col_indexes].tolist()
        return pd.DataFrame(
            [src_obj_ids, target_obj_ids, scores_candidate], index=["source", "target", "score"]
        )

    def _calc_standardised_score_met(self, results: list) -> list[pd.DataFrame]:
        if self.metcalf_mean is None or self.metcalf_std is None:
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
            mean = self.metcalf_mean[num_met_strains][num_gcf_strains]
            sqrt = self.metcalf_std[num_met_strains][num_gcf_strains]
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

    def _calc_standardised_score_gen(self, results: list) -> list[pd.DataFrame]:
        if self.metcalf_mean is None or self.metcalf_std is None:
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
                mean = self.metcalf_mean[num_met_strains][num_gcf_strains]
                sqrt = self.metcalf_std[num_met_strains][num_gcf_strains]
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
