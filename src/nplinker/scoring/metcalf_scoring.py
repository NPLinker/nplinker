from __future__ import annotations
import logging
from enum import Enum
from typing import TYPE_CHECKING
from typing import TypeVar
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from nplinker.genomics import GCF
from nplinker.metabolomics import MolecularFamily
from nplinker.metabolomics import Spectrum
from .abc import ScoringBase
from .link_graph import LinkGraph
from .link_graph import Score
from .utils import get_presence_gcf_strain
from .utils import get_presence_mf_strain
from .utils import get_presence_spec_strain


if TYPE_CHECKING:
    from ..nplinker import NPLinker


logger = logging.getLogger(__name__)


class LinkType(Enum):
    """Enum class for link types."""

    SPEC_GCF = "spec-gcf"
    MF_GCF = "mf-gcf"


ObjectType = TypeVar("ObjectType", GCF, Spectrum, MolecularFamily)


class MetcalfScoring(ScoringBase):
    """Metcalf scoring method.

    Attributes:
        name: The name of this scoring method, set to a fixed value `metcalf`.
        npl: The NPLinker object.
        CACHE: The name of the cache file to use for storing the MetcalfScoring.

        presence_gcf_strain: A DataFrame to store presence of gcfs with respect to strains.
            The index of the DataFrame are the GCF objects and the columns are Strain objects.
            The values are 1 where the gcf occurs in the strain, 0 otherwise.
        presence_spec_strain: A DataFrame to store presence of spectra with respect to strains.
            The index of the DataFrame are the Spectrum objects and the columns are Strain objects.
            The values are 1 where the spectrum occurs in the strain, 0 otherwise.
        presence_mf_strain: A DataFrame to store presence of molecular families with respect to strains.
            The index of the DataFrame are the MolecularFamily objects and the columns are Strain objects.
            The values are 1 where the molecular family occurs in the strain, 0 otherwise.

        raw_score_spec_gcf: A DataFrame to store the raw Metcalf scores for spectrum-gcf links.
            The columns are "spec", "gcf" and "score":

            - The "spec" and "gcf" columns contain the Spectrum and GCF objects respectively,
            - The "score" column contains the raw Metcalf scores.

        raw_score_mf_gcf: A DataFrame to store the raw Metcalf scores for molecular family-gcf links.
            The columns are "mf", "gcf" and "score":

            - The "mf" and "gcf" columns contain the MolecularFamily and GCF objects respectively,
            - the "score" column contains the raw Metcalf scores.

        metcalf_mean: The mean value used for standardising Metcalf scores.
        metcalf_std: The standard deviation value used for standardising Metcalf scores.
    """

    name = "metcalf"
    npl: NPLinker | None = None
    CACHE: str = "cache_metcalf_scoring.pckl"
    metcalf_weights: tuple[int, int, int, int] = (10, -10, 0, 1)

    # index: gcf/spec/mf ids, columns: strain ids, value: 0/1
    presence_gcf_strain: pd.DataFrame = pd.DataFrame()
    presence_spec_strain: pd.DataFrame = pd.DataFrame()
    presence_mf_strain: pd.DataFrame = pd.DataFrame()

    raw_score_spec_gcf: pd.DataFrame = pd.DataFrame(columns=["spec", "gcf", "score"])
    raw_score_mf_gcf: pd.DataFrame = pd.DataFrame(columns=["mf", "gcf", "score"])

    metcalf_mean: np.ndarray | None = None
    metcalf_std: np.ndarray | None = None

    @classmethod
    def setup(cls, npl: NPLinker):
        """Setup the MetcalfScoring object.

        This method is only called once to setup the MetcalfScoring object.

        Args:
            npl: The NPLinker object.
        """
        if cls.npl is not None:
            logger.info("MetcalfScoring.setup already called, skipping.")
            return

        logger.info(
            f"MetcalfScoring.setup starts: #bgcs={len(npl.bgcs)}, #gcfs={len(npl.gcfs)}, "
            f"#spectra={len(npl.spectra)}, #molfams={len(npl.molfams)}, #strains={npl.strains}"
        )
        cls.npl = npl

        # calculate presence of gcfs/spectra/mfs with respect to strains
        cls.presence_gcf_strain = get_presence_gcf_strain(npl.gcfs, npl.strains)
        cls.presence_spec_strain = get_presence_spec_strain(npl.spectra, npl.strains)
        cls.presence_mf_strain = get_presence_mf_strain(npl.molfams, npl.strains)

        # calculate raw Metcalf scores for spec-gcf links
        raw_score_spec_gcf = cls._calc_raw_score(
            cls.presence_spec_strain, cls.presence_gcf_strain, cls.metcalf_weights
        )
        cls.raw_score_spec_gcf = raw_score_spec_gcf.reset_index().melt(id_vars="index")
        cls.raw_score_spec_gcf.columns = ["spec", "gcf", "score"]

        # calculate raw Metcalf scores for spec-gcf links
        raw_score_mf_gcf = cls._calc_raw_score(
            cls.presence_mf_strain, cls.presence_gcf_strain, cls.metcalf_weights
        )
        cls.raw_score_mf_gcf = raw_score_mf_gcf.reset_index().melt(id_vars="index")
        cls.raw_score_mf_gcf.columns = ["mf", "gcf", "score"]

        # calculate mean and std for standardising Metcalf scores
        n_strains = cls.presence_gcf_strain.shape[1]
        cls.metcalf_mean, cls.metcalf_std = cls._calc_mean_std(n_strains, cls.metcalf_weights)

        logger.info("MetcalfScoring.setup completed")

    def get_links(self, *objects: ObjectType, **parameters) -> LinkGraph:
        """Get links for the given objects.

        Args:
            objects: The objects to get links for. All objects must be of the same type, i.e. `GCF`,
              `Spectrum` or `MolecularFamily` type.
               If no objects are provided, all detected objects (`npl.gcfs`) will be used.
            parameters: The scoring parameters to use for the links. The parameters are:

                    - cutoff: The minimum score to consider a link (≥cutoff). Default is 0.
                    - standardised: Whether to use standardised scores. Default is False.

        Returns:
            The `LinkGraph` object containing the links involving the input objects with the Metcalf
            scores.

        Raises:
            TypeError: If the input objects are not of the same type or the object type is invalid.
        """
        # validate input objects
        if len(objects) == 0:
            objects = self.npl.gcfs
        # check if all objects are of the same type
        types = {type(i) for i in objects}
        if len(types) > 1:
            raise TypeError("Input objects must be of the same type.")
        # check if the object type is valid
        obj_type = next(iter(types))
        if obj_type not in (GCF, Spectrum, MolecularFamily):
            raise TypeError(
                f"Invalid type {obj_type}. Input objects must be GCF, Spectrum or MolecularFamily objects."
            )

        # validate scoring parameters
        self._cutoff: float = parameters.get("cutoff", 0)
        self._standardised: bool = parameters.get("standardised", False)
        parameters.update({"cutoff": self._cutoff, "standardised": self._standardised})

        logger.info(
            f"MetcalfScoring: #objects={len(objects)}, type={obj_type}, cutoff={self._cutoff}, "
            f"standardised={self._standardised}"
        )
        if not self._standardised:
            scores_list = self._get_links(*objects, obj_type=obj_type, score_cutoff=self._cutoff)
        else:
            # use negative infinity as the score cutoff to ensure we get all links
            scores_list = self._get_links(*objects, obj_type=obj_type, score_cutoff=np.NINF)
            if obj_type == GCF:
                scores_list = self._calc_standardised_score_gen(scores_list)
            else:
                scores_list = self._calc_standardised_score_met(scores_list)

        links = LinkGraph()
        for score_df in scores_list:
            for row in score_df.itertuples(index=False):  # row has attributes: spec/mf, gcf, score
                met = row.spec if score_df.name == LinkType.SPEC_GCF else row.mf
                links.add_link(
                    row.gcf,
                    met,
                    metcalf=Score(self.name, row.score, parameters),
                )

        logger.info(f"MetcalfScoring: completed! Found {len(links.links)} links in total.")
        return links

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
    def _calc_raw_score(
        p1: pd.DataFrame, p2: pd.DataFrame, weights: tuple[int, int, int, int]
    ) -> pd.DataFrame:
        """Calculate non-standardised Metcalf scores.

        Args:
            p1: A DataFrame containing the presence of objects in strains.
            p2: A DataFrame containing the presence of objects in strains.
            weights: The weights to use for Metcalf scoring.

        Returns:
            A DataFrame containing the non-standardised Metcalf scores.
        """
        nop1 = 1 - p1
        nop2 = 1 - p2

        # calculate co-presence
        p1_p2 = p1.dot(p2.T)
        p1_nop2 = p1.dot(nop2.T)
        nop1_p2 = nop1.dot(p2.T)
        nop1_nop2 = nop1.dot(nop2.T)

        # calculate weighted sum
        score = (
            p1_p2 * weights[0]
            + p1_nop2 * weights[1]
            + nop1_p2 * weights[2]
            + nop1_nop2 * weights[3]
        )

        return score

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
        *objects: ObjectType,
        obj_type: GCF | Spectrum | MolecularFamily,
        score_cutoff: float = 0,
    ) -> list[pd.DataFrame]:
        """Get links and scores for the given objects.

        Args:
            objects: A list of GCF, Spectrum or MolecularFamily objects and all objects must be of
                the same type.
            obj_type: The type of the objects.
            score_cutoff: Minimum score to consider a link (≥score_cutoff). Default is 0.

        Returns:
            List of data frames containing the ids of the linked objects and the score.

            The data frame is named by link types, see `LinkType`. It has column names of
            ['spec', 'gcf', 'score'] or ['mf', 'gcf', 'score'] depending on the link type:

            - the 'spec', 'mf' or 'gcf' column contains the Spectrum, MolecularFamily or GCF objects,
            - the 'score' column contains the scores.
        """
        links = []
        col_name = (
            "spec" if obj_type == Spectrum else "mf" if obj_type == MolecularFamily else "gcf"
        )

        # spec-gcf link
        if obj_type in (GCF, Spectrum):
            df = self.raw_score_spec_gcf[
                self.raw_score_spec_gcf[col_name].isin(objects)
                & (self.raw_score_spec_gcf["score"] >= score_cutoff)
            ]
            df.name = LinkType.SPEC_GCF
            links.append(df)

        # mf-gcf link
        if obj_type in (GCF, MolecularFamily):
            df = self.raw_score_mf_gcf[
                self.raw_score_mf_gcf[col_name].isin(objects)
                & (self.raw_score_mf_gcf["score"] >= score_cutoff)
            ]
            df.name = LinkType.MF_GCF
            links.append(df)

        return links

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
            if raw_score.name == LinkType.SPEC_GCF:
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
        mask = z_scores >= self._cutoff

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
                if raw_score.name == LinkType.SPEC_GCF:
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
            mask = z_scores >= self._cutoff

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
