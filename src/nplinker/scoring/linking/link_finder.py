from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from nplinker.genomics.gcf import GCF
from nplinker.logconfig import LogConfig
from nplinker.metabolomics import MolecularFamily
from nplinker.metabolomics import Spectrum
from . import LINK_TYPES
from .utils import isinstance_all


if TYPE_CHECKING:
    from . import DataLinks

logger = LogConfig.getLogger(__file__)


# TODO CG: this class could be merged to MetcalfScoring class?
class LinkFinder:
    def __init__(self) -> None:
        """Initialise LinkFinder object.

        Attributes:
            raw_score_spec_gcf: The raw Metcalf scores for
                spectrum-GCF links.
            raw_score_mf_gcf: The raw Metcalf scores for
                molecular family-GCF links.
            metcalf_mean: The mean value used for
                standardising Metcalf scores.
            metcalf_std: The standard deviation value used
                for standardising Metcalf scores.
        """
        self.raw_score_spec_gcf = pd.DataFrame()
        self.raw_score_mf_gcf = pd.DataFrame()
        self.metcalf_mean = None
        self.metcalf_std = None

    # TODO CG: calc_score method could be integrated to __init__?
    def calc_score(
        self,
        data_links: DataLinks,
        link_type: str = "spec-gcf",
        scoring_weights: tuple[int, int, int, int] = (10, -10, 0, 1),
    ) -> None:
        """Calculate Metcalf scores.

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
            self.raw_score_spec_gcf = (
                data_links.cooccurrence_spec_gcf * scoring_weights[0]
                + data_links.cooccurrence_spec_notgcf * scoring_weights[1]
                + data_links.cooccurrence_notspec_gcf * scoring_weights[2]
                + data_links.cooccurrence_notspec_notgcf * scoring_weights[3]
            )
        if link_type == "mf-gcf":
            self.raw_score_mf_gcf = (
                data_links.cooccurrence_mf_gcf * scoring_weights[0]
                + data_links.cooccurrence_mf_notgcf * scoring_weights[1]
                + data_links.cooccurrence_notmf_gcf * scoring_weights[2]
                + data_links.cooccurrence_notmf_notgcf * scoring_weights[3]
            )

        # TODO CG: this part should be moved outside of this method
        n_strains = data_links.occurrence_gcf_strain.shape[1]
        if self.metcalf_mean is None or self.metcalf_std is None:
            self.metcalf_mean, self.metcalf_std = self._calc_mean_std(n_strains, scoring_weights)

    # TODO CG: read paper and check the logics of this method
    def _calc_mean_std(
        self, n_strains: int, scoring_weights: tuple[int, int, int, int]
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

    def get_links(
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

    def _get_scores_source_gcf(self, scores: pd.DataFrame, score_cutoff: float) -> pd.DataFrame:
        row_indexes, col_indexes = np.where(scores >= score_cutoff)
        src_obj_ids = scores.columns[col_indexes].to_list()
        target_obj_ids = scores.index[row_indexes].to_list()
        scores_candidate = scores.values[row_indexes, col_indexes].tolist()
        return pd.DataFrame(
            [src_obj_ids, target_obj_ids, scores_candidate], index=["source", "target", "score"]
        )

    def _get_scores_source_met(self, scores: pd.DataFrame, score_cutoff: float) -> pd.DataFrame:
        row_indexes, col_indexes = np.where(scores >= score_cutoff)
        src_obj_ids = scores.index[row_indexes].to_list()
        target_obj_ids = scores.columns[col_indexes].to_list()
        scores_candidate = scores.values[row_indexes, col_indexes].tolist()
        return pd.DataFrame(
            [src_obj_ids, target_obj_ids, scores_candidate], index=["source", "target", "score"]
        )
