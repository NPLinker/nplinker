from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from nplinker.genomics.gcf import GCF
from nplinker.logconfig import LogConfig
from nplinker.metabolomics.molecular_family import MolecularFamily
from nplinker.metabolomics.spectrum import Spectrum


if TYPE_CHECKING:
    from .data_linking import DataLinks

logger = LogConfig.getLogger(__file__)

LINK_TYPES = ['spec-gcf', 'mf-gcf']


class LinkFinder():
    """
    Class to:
    1) Score potential links based on collected information from:
        DataLinks, LinkLikelihood (and potentially other resources)
    Different scores can be used for this!

    2) Rank and output selected candidates
    3) Create output plots and tables
    """

    def __init__(self):
        """
        Create tables of prospective link candidates.
        Separate tables will exist for different linking scenarios, such as
        gcfs <-> spectra OR gcf <-> mol.families
        """
        # metcalf scores
        self.raw_score_spec_gcf = pd.DataFrame()
        self.raw_score_fam_gcf = pd.DataFrame()

        # metcalf caching
        self.metcalf_mean = None
        self.metcalf_std = None

    def cal_score(
        self,
        data_links: DataLinks,
        link_type: str = 'spec-gcf',
        scoring_weights: tuple[int, int, int, int] = (10, -10, 0, 1)
    ) -> None:
        """Calculate metcalf scores.

        Args:
            data_links (DataLinks): The DataLinks object to use for scoring.
            link_type (str, optional): The type of link to score. Available
                types are 'spec-gcf' and 'mf-gcf'. Defaults to 'spec-gcf'.
            scoring_weights (tuple[int,int,int,int], optional): The weights to
                use for Metcalf scoring. The weights are applied to
                '(met_gcf, met_not_gcf, gcf_not_met, not_met_not_gcf)'.
                Defaults to (10, -10, 0, 1).
        """
        if link_type == 'spec-gcf':
            self.raw_score_spec_gcf = (
                data_links.cooccurrence_spec_gcf * scoring_weights[0] +
                data_links.cooccurrence_spec_notgcf * scoring_weights[1] +
                data_links.cooccurrence_notspec_gcf * scoring_weights[2] +
                data_links.cooccurrence_notspec_notgcf * scoring_weights[3])
        if link_type == 'mf-gcf':
            self.raw_score_fam_gcf = (
                data_links.cooccurrence_mf_gcf * scoring_weights[0] +
                data_links.cooccurrence_mf_notgcf * scoring_weights[1] +
                data_links.cooccurrence_notmf_gcf * scoring_weights[2] +
                data_links.cooccurrence_notmf_notgcf * scoring_weights[3])

        if self.metcalf_mean is None or self.metcalf_std is None:
            self.metcalf_mean, self.metcalf_std = self._cal_mean_std(
                data_links, scoring_weights)

    def _cal_mean_std(self, data_links, scoring_weights):
        # Compute the expected values for all possible values of spec and gcf strains
        # we need the total number of strains
        _, n_strains = data_links.occurrence_gcf_strain.shape
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
        *objects: tuple[GCF, ...] | tuple[Spectrum, ...]
        | tuple[MolecularFamily, ...],
        score_cutoff: float = 0.5
    ) -> list[tuple[tuple[str], tuple[str], tuple[float]]]:
        """Get scores for links between objects.

        Args:
            objects(tuple): GCF, Spectrum or MolecularFamily objects.
            score_cutoff(float): Minimum score to consider a link (≥score_cutoff).
                Default is 0.5.
        Returns:
            list: List of tuples containing the ids of the linked objects and the score.
                The tuple contains three tuples:
                - the first tuple contains the ids of the input/source objects,
                - the second tuple contains the ids of the target objects,
                - the third tuple contains the scores.

        Raises:
            TypeError: If input objects are not GCF, Spectrum or MolecularFamily objects.
        """
        if self._isinstance(*objects, GCF):
            obj_type = 'gcf'
        elif self._isinstance(*objects, Spectrum):
            obj_type = 'spec'
        elif self._isinstance(*objects, MolecularFamily):
            obj_type = 'fam'
        else:
            raise TypeError(
                'Input objects must be GCF, Spectrum or MolecularFamily objects.'
            )

        links = []
        if obj_type == 'gcf':
            obj_ids = [gcf.gcf_id for gcf in objects]
            all_scores = (self.raw_score_spec_gcf.loc[:, obj_ids],
                          self.raw_score_fam_gcf.loc[:, obj_ids])
            for scores in all_scores:
                links.append(self._get_scores_source_gcf(scores, score_cutoff))
        if obj_type == 'spec':
            obj_ids = [spec.spectrum_id for spec in objects]
            all_scores = self.raw_score_spec_gcf.loc[obj_ids, :]
            links.append(self._get_scores_source_met(all_scores, score_cutoff))
        if obj_type == 'fam':
            obj_ids = [mf.family_id for mf in objects]
            all_scores = self.raw_score_fam_gcf.loc[obj_ids, :]
            links.append(self._get_scores_source_met(all_scores, score_cutoff))
        return links

    def _isinstance(self, *objects, obj_type=GCF) -> bool:
        return all(isinstance(x, obj_type) for x in objects)

    def _get_scores_source_gcf(self, scores, score_cutoff):
        candidate_met_gcf_indexes = np.where(scores >= score_cutoff)
        src_obj_ids = scores.columns[candidate_met_gcf_indexes[1]].to_list()
        target_obj_ids = scores.index[candidate_met_gcf_indexes[0]].to_list()
        scores_candidate = scores.to_numpy()[candidate_met_gcf_indexes].tolist(
        )
        return tuple(src_obj_ids), tuple(target_obj_ids), tuple(
            scores_candidate)

    def _get_scores_source_met(self, scores, score_cutoff):
        candidate_met_gcf_indexes = np.where(scores >= score_cutoff)
        src_obj_ids = scores.index[candidate_met_gcf_indexes[0]].to_list()
        target_obj_ids = scores.columns[candidate_met_gcf_indexes[1]].to_list()
        scores_candidate = scores.to_numpy()[candidate_met_gcf_indexes].tolist(
        )
        return tuple(src_obj_ids), tuple(target_obj_ids), tuple(
            scores_candidate)
