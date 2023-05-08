from __future__ import annotations
from typing import TYPE_CHECKING
from deprecated import deprecated
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from nplinker.genomics.gcf import GCF
from nplinker.metabolomics.spectrum import Spectrum
from .data_linking_functions import pair_prob_approx
from .data_linking_functions import pair_prob_hg


# CG: TODO get_links function does not work any more, need to update its logics

# import packages for plotting
# TODO move plotting to separate module?
try:
    from matplotlib import pyplot as plt
    import seaborn as sns
except ImportError:
    print(
        'Warning: plotting functionality will not be available (missing matplotlib and/or seaborn)'
    )

from nplinker.logconfig import LogConfig
from nplinker.metabolomics.molecular_family import MolecularFamily


if TYPE_CHECKING:
    from .data_linking import DataLinks

logger = LogConfig.getLogger(__file__)

SCORE_TYPES = ['metcalf', 'likescore', 'hg']
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
        self.metcalf_spec_gcf = pd.DataFrame()
        self.metcalf_fam_gcf = pd.DataFrame()

        # likelihood scores
        self.likescores_spec_gcf = []
        self.likescores_fam_gcf = []

        # hg scores
        self.hg_spec_gcf = []
        self.hg_fam_gcf = []

        # link candidate tables
        self.link_candidates_gcf_spec = []
        self.link_candidates_gcf_fam = []

        # metcalf caching
        self.metcalf_mean = None
        self.metcalf_std = None

    def get_scores(self, score_type: str, link_type: str) -> pd.DataFrame:
        """Get the scores for given scoring type and link type.

        Args:
            score_type (str): The type of scoring method to use. Available
                scoring methods are 'metcalf', 'likescore', and 'hg'.
            link_type (str): The type of link to get scores for. Available
                types are 'spec-gcf' and 'mf-gcf'.

        Returns:
            pd.DataFrame: The scores for the given method and link type.

        Raises:
            ValueError: If an invalid score type or link type is given.
        """
        if score_type not in SCORE_TYPES:
            raise ValueError(
                f'Invalid score type "{score_type}". Must be one of "{SCORE_TYPES}"'
            )
        if link_type not in LINK_TYPES:
            raise ValueError(
                f'Invalid link type "{link_type}". Must be one of "{LINK_TYPES}"'
            )
        if score_type == 'metcalf':
            if link_type == 'spec-gcf':
                return self.metcalf_spec_gcf
            if link_type == 'mf-gcf':
                return self.metcalf_fam_gcf
        if score_type == 'likescore':
            if link_type == 'spec-gcf':
                return self.likescores_spec_gcf
            if link_type == 'mf-gcf':
                return self.likescores_fam_gcf
        if score_type == 'hg':
            if link_type == 'spec-gcf':
                return self.hg_spec_gcf
            if link_type == 'mf-gcf':
                return self.hg_fam_gcf

    def metcalf_scoring(
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
            self.metcalf_spec_gcf = (
                data_links.cooccurrence_spec_gcf * scoring_weights[0] +
                data_links.cooccurrence_spec_notgcf * scoring_weights[1] +
                data_links.cooccurrence_notspec_gcf * scoring_weights[2] +
                data_links.cooccurrence_notspec_notgcf * scoring_weights[3])
        if link_type == 'mf-gcf':
            self.metcalf_fam_gcf = (
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

    def hg_scoring(self, data_links, type='spec-gcf'):
        """
        Calculate metcalf scores from DataLinks() co-occurence matrices
        """

        # NOTE:can't use the correlation matrices directly for this scoring method because
        # it seems to require more inclusive counts of the strains in each object.
        # Instead of "number of strains only in GCF", it requires "number of strains in the
        # GCF PLUS the number shared between the GCF and the other object".
        # e.g. if a spectrum has 3 strains, a GCF has 1 strain and there is 1 shared strain,
        # cooccurrence_spec_gcf will correctly contain "1", but M_type1_notgcf will contain "2" instead
        # of "3", because the spectrum only has 2 distinct strains vs the GCF.
        # To fix this the cooccurrence_spec_gcf/cooccurrence_mf_gcf matrix can just be added onto the others to give
        # the correct totals.

        if type == 'spec-gcf':
            num_strains = np.ones(data_links.cooccurrence_spec_gcf.shape
                                  ) * data_links.occurrence_gcf_strain.shape[1]
            overlap_counts = data_links.cooccurrence_spec_gcf
            gcf_counts = overlap_counts + data_links.cooccurrence_notspec_gcf
            spec_counts = overlap_counts + data_links.cooccurrence_spec_notgcf
            hg_scores = hypergeom.sf(overlap_counts,
                                     num_strains,
                                     gcf_counts,
                                     spec_counts,
                                     loc=1)
            self.hg_spec_gcf = hg_scores
        elif type == 'mf-gcf':
            num_strains = np.ones(data_links.cooccurrence_mf_gcf.shape
                                  ) * data_links.occurrence_gcf_strain.shape[1]
            overlap_counts = data_links.cooccurrence_mf_gcf
            gcf_counts = overlap_counts + data_links.cooccurrence_notmf_gcf
            fam_counts = overlap_counts + data_links.cooccurrence_mf_notgcf
            hg_scores = hypergeom.sf(overlap_counts,
                                     num_strains,
                                     gcf_counts,
                                     fam_counts,
                                     loc=1)
            self.hg_fam_gcf = hg_scores

        return hg_scores

    def likelihood_scoring(self,
                           data_links,
                           likelihoods,
                           alpha_weighing=0.5,
                           type='spec-gcf'):
        """
        Calculate likelihood scores from DataLinks() co-occurence matrices.

        Idea:
            Score reflect the directionality BGC-->compound-->spectrum, which
            suggests that the most relevant likelihoods are:
            P(gcf|type1) - If type1 is the result of only one particular gene cluster,
                            this value should be high (close or equal to 1)
            P(type1|not gcf) - Following the same logic, this value should be very
                                small or 0 ("no gene cluster, no compound")

        Score:
            Score = P(gcf|type1) * (1 - P(type1|not gcf) * weighing function

            weighing function is here a function of the number of strains they co-occur.
            weighing function = (1 - exp(-alpha_weighing * num_of_co-occurrences)
        """

        if type == 'spec-gcf':
            likelihood_scores = np.zeros(
                data_links.cooccurrence_spec_gcf.shape)
            likelihood_scores = (
                likelihoods.P_gcf_given_spec *
                (1 - likelihoods.P_spec_not_gcf) *
                (1 -
                 np.exp(-alpha_weighing * data_links.cooccurrence_spec_gcf)))

            self.likescores_spec_gcf = likelihood_scores

        elif type == 'mf-gcf':
            likelihood_scores = np.zeros(data_links.cooccurrence_mf_gcf.shape)
            likelihood_scores = (
                likelihoods.P_gcf_given_fam * (1 - likelihoods.P_fam_not_gcf) *
                (1 - np.exp(-alpha_weighing * data_links.cooccurrence_mf_gcf)))

            self.likescores_fam_gcf = likelihood_scores
        return likelihood_scores

    def select_link_candidates(self,
                               data_links,
                               likelihoods,
                               P_cutoff=0.8,
                               main_score='likescore',
                               score_cutoff=0,
                               type='mf-gcf'):
        """
        Look for potential best candidate for links between
        IF type='spec-gcf': GCFs and spectra
        IF type='mf-gcf': GCFs and mol.families

        Parameters
        ----------
        data_links: DataLinks() object
            containing co-occurence matrices
        likelihood: LinkLikelihood() object
            containing co-occurence likelihoods
        P_cutoff: float
            Thresholds to conly consider candidates for which:
            P_gcf_given_type1 >= P_cutoff
        main_score: str
            Which main score to use ('metcalf', 'likescore')
        score_cutoff:
            Thresholds to conly consider candidates for which:
            score >= score_cutoff
        """

        # Select scenario: spec<->gcf or mf<->gcf
        if type == 'spec-gcf':
            P_gcf_given_type1 = likelihoods.P_gcf_given_spec
            P_gcf_not_type1 = likelihoods.P_gcf_not_spec
            P_type1_given_gcf = likelihoods.P_spec_given_gcf
            P_type1_not_gcf = likelihoods.P_spec_not_gcf
            M_type1_gcf = data_links.cooccurrence_spec_gcf
            metcalf_scores = self.metcalf_spec_gcf
            likescores = self.likescores_spec_gcf
            index_names = [
                "spectrum_id", "GCF id", "P(gcf|spec)", "P(spec|gcf)",
                "P(gcf|not spec)", "P(spec|not gcf)", "co-occur in # strains",
                "metcalf score", "likelihood score", "HG prob", "link prob",
                "link prob specific"
            ]

        elif type == 'mf-gcf':
            P_gcf_given_type1 = likelihoods.P_gcf_given_fam
            P_gcf_not_type1 = likelihoods.P_gcf_not_fam
            P_type1_given_gcf = likelihoods.P_fam_given_gcf
            P_type1_not_gcf = likelihoods.P_fam_not_gcf
            M_type1_gcf = data_links.cooccurrence_mf_gcf
            metcalf_scores = self.metcalf_fam_gcf
            likescores = self.likescores_fam_gcf
            index_names = [
                "family_id", "GCF id", "P(gcf|fam)", "P(fam|gcf)",
                "P(gcf|not fam)", "P(fam|not gcf)", "co-occur in # strains",
                "metcalf score", "likelihood score", "HG prob", "link prob",
                "link prob specific"
            ]

        elif type == 'spec-bgc' or type == 'fam-bgc':
            raise Exception("Given types are not yet supported... ")
        else:
            raise Exception(
                "Wrong correlation 'type' given. Must be one of 'spec-gcf', 'mf-gcf'..."
            )

        dim1, dim2 = P_gcf_given_type1.shape

        # PRE-SELECTION:
        # Select candidates with P_gcf_given_spec >= P_cutoff AND score >= score_cutoff
        if main_score == 'likescore':
            candidate_ids = np.where((P_gcf_given_type1[:, :] >= P_cutoff)
                                     & (likescores >= score_cutoff))
        elif main_score == 'metcalf':
            candidate_ids = np.where((P_gcf_given_type1[:, :] >= P_cutoff)
                                     & (metcalf_scores >= score_cutoff))
        else:
            raise Exception(
                'Wrong scoring type given. Must be one of: {}'.format(
                    SCORE_TYPES))

        logger.debug(candidate_ids[0].shape[0], " candidates selected with ",
                     index_names[2], " >= ", P_cutoff, " and a link score >= ",
                     score_cutoff, ".")

        link_candidates = np.zeros((12, candidate_ids[0].shape[0]))
        link_candidates[0, :] = candidate_ids[0]  # spectrum/fam number
        link_candidates[1, :] = candidate_ids[1]  # gcf id
        link_candidates[2, :] = P_gcf_given_type1[candidate_ids]
        link_candidates[3, :] = P_type1_given_gcf[candidate_ids]
        link_candidates[4, :] = P_gcf_not_type1[candidate_ids]
        link_candidates[5, :] = P_type1_not_gcf[candidate_ids]
        link_candidates[6, :] = M_type1_gcf[candidate_ids]
        link_candidates[7, :] = metcalf_scores[candidate_ids]
        link_candidates[8, :] = likescores[candidate_ids]

        # Calculate probability to find similar link by chance
        Nx_list = data_links.mapping_gcf["no of strains"]
        if type == 'spec-gcf':
            Ny_list = data_links.mapping_spec["no of strains"]
        elif type == 'mf-gcf':
            Ny_list = data_links.mapping_fam["no of strains"]

        # Calculate probabilities of finding a spectrum in a certain strain
        P_str = np.array(data_links.mapping_strain["no of spectra"])
        P_str = P_str / np.sum(P_str)

        num_strains = data_links.occurrence_gcf_strain.shape[1]

        # Calculate the hypergeometric probability (as before)
        for i in range(link_candidates.shape[1]):
            link_candidates[9, i] = pair_prob_hg(
                link_candidates[6, i], num_strains,
                Nx_list[link_candidates[1, i]],
                Ny_list[int(link_candidates[0, i])])

        # Calculate the GCF specific probability
        for i in range(link_candidates.shape[1]):
            id_spec = link_candidates[0, i]
            id_gcf = link_candidates[1, i]

            # find set of strains which contain GCF with id link_candidates[1,i]
            XG = np.where(
                data_links.occurrence_gcf_strain.loc[id_gcf, :] == 1)[0]

            link_candidates[10,
                            i] = pair_prob_approx(P_str, XG,
                                                  int(Ny_list[id_spec]),
                                                  int(link_candidates[6, i]))
            # Calculate the link specific probability
            # Find strains where GCF and spectra/family co-occur
            if type == 'spec-gcf':
                XGS = np.where(
                    (data_links.occurrence_gcf_strain[id_gcf, :] == 1)
                    & (data_links.occurrence_spec_strain[id_spec, :] == 1))[0]
            elif type == 'mf-gcf':
                XGS = np.where(
                    (data_links.occurrence_gcf_strain[id_gcf, :] == 1)
                    & (data_links.occurrence_mf_strain[id_spec, :] == 1))[0]
            link_candidates[11,
                            i] = link_prob(P_str, XGS, int(Nx_list[id_gcf]),
                                           int(Ny_list[id_spec]), num_strains)

        # Transform into pandas Dataframe (to avoid index confusions):
        link_candidates_pd = pd.DataFrame(link_candidates.transpose(1, 0),
                                          columns=index_names)

        # add other potentially relevant knowdledge
        # If this will grow to more collected information -> create separate function/class
        bgc_class = []
        # TODO CG: bgc class should be obtained from GCF object
        for i in link_candidates_pd["GCF id"].astype(int):
            bgc_class.append(data_links.mapping_gcf["bgc class"][i])
        link_candidates_pd["BGC class"] = bgc_class

        # Change some columns to int
        link_candidates_pd.iloc[:, [0, 1, 6, 7]] = link_candidates_pd.iloc[:, [
            0, 1, 6, 7
        ]].astype(int)

        # return results
        if type == 'spec-gcf':
            self.link_candidates_gcf_spec = link_candidates_pd
        elif type == 'mf-gcf':
            self.link_candidates_gcf_fam = link_candidates_pd
        else:
            raise Exception("No candidate selection was created.")
        return link_candidates_pd

    @staticmethod
    def _isinstance(*objects, obj_type=GCF) -> bool:
        return all(isinstance(x, obj_type) for x in objects)

    def get_links(
        self,
        *objects: tuple[GCF, ...] | tuple[Spectrum, ...]
        | tuple[MolecularFamily, ...],
        score_type: str = 'likescore',
        score_cutoff: float = 0.5
    ) -> list[tuple[tuple[str], tuple[str], tuple[float]]]:
        """Get scores for links between objects.

        Args:
            objects(tuple): GCF, Spectrum or MolecularFamily objects.
            score_type(str): Type of score to use for link calculation.
                Must be one of 'likescore', 'metcalf' and 'hg'.
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
            ValueError: If score type is not one of 'likescore', 'metcalf' and 'hg'.
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

        if score_type not in SCORE_TYPES:
            raise ValueError(
                f'Invalid score type "{score_type}". Must be one of "{SCORE_TYPES}"'
            )

        links = []
        # TODO CG: replace integer ids with string ids for `hg` and `likescore`
        if obj_type == 'gcf':
            obj_ids = [gcf.gcf_id for gcf in objects]
            if score_type == 'likescore':
                all_scores = (self.likescores_spec_gcf[:, obj_ids],
                              self.likescores_fam_gcf[:, obj_ids])
            if score_type == 'metcalf':
                all_scores = (self.metcalf_spec_gcf.loc[:, obj_ids],
                              self.metcalf_fam_gcf.loc[:, obj_ids])
            if score_type == 'hg':
                all_scores = (self.hg_spec_gcf[:, obj_ids],
                              self.hg_fam_gcf[:, obj_ids])
            for scores in all_scores:
                links.append(self._get_scores_source_gcf(scores, score_cutoff))

        if obj_type == 'spec':
            obj_ids = [spec.spectrum_id for spec in objects]
            if score_type == 'likescore':
                all_scores = self.likescores_spec_gcf[obj_ids, :]
            if score_type == 'metcalf':
                all_scores = self.metcalf_spec_gcf.loc[obj_ids, :]
            if score_type == 'hg':
                all_scores = self.hg_spec_gcf[obj_ids, :]
            links.append(self._get_scores_source_met(all_scores, score_cutoff))

        if obj_type == 'fam':
            obj_ids = [mf.family_id for mf in objects]
            if score_type == 'likescore':
                all_scores = self.likescores_fam_gcf[obj_ids, :]
            if score_type == 'metcalf':
                all_scores = self.metcalf_fam_gcf.loc[obj_ids, :]
            if score_type == 'hg':
                all_scores = self.hg_fam_gcf[obj_ids:]
            links.append(self._get_scores_source_met(all_scores, score_cutoff))

        return links

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

    @deprecated(version="1.3.3",
                reason="The unworkable method will be removed")
    def create_cytoscape_files(self,
                               data_links,
                               network_filename,
                               link_type='mf-gcf',
                               score_type='metcalf'):
        """
        Create network file for import into Cytoscape.
        Network file will be generated using networkx.
        The type of network created here is a bipartite network.
            mass spec side --> bipartite = 0
            gene cluster side --> bipartite = 1
        Output format is a graphml file.

        Parameters
        ----------
        data_links: DataLinks() object
            containing co-occurence matrices
        likelihood: LinkLikelihood() object
            containing co-occurence likelihoods
        network_filename: str
            Filename to save generated model as graphml file.
        """

        import networkx as nx
        NPlinker_net = nx.Graph()

        if link_type == 'mf-gcf':
            link_candidates = self.link_candidates_gcf_fam
            type1str = 'family_id'
        elif link_type == 'spec-gcf':
            link_candidates = self.link_candidates_gcf_spec
            type1str = 'spectrum_id'
        else:
            raise Exception("Wrong link-type given.")

        # Select score type
        if score_type == 'metcalf':
            scorestr = 'metcalf score'
        elif score_type == 'likescore':
            scorestr = 'likelihood score'
        else:
            raise Exception(
                "Wrong score_type given. Must be one of: 'metcalf', 'likescore' ."
            )

        # Add nodes (all partners from link_candidate table):
        # mass spec side --> bipartite = 0
        type1_names = []
        for type1 in link_candidates[type1str].astype(int):
            type1_names.append(type1str + str(type1))

        type1_names_unique = list(set(type1_names))
        NPlinker_net.add_nodes_from(type1_names_unique, bipartite=0)

        # gene cluster side --> bipartite = 1
        type2_names = []
        for type2 in link_candidates['GCF id']:
            type2_names.append("GCF_" + str(type2))

        type2_names_unique = list(set(type2_names))
        NPlinker_net.add_nodes_from(type2_names_unique, bipartite=1)

        # Add edges:
        for i in range(0, link_candidates.shape[0]):
            NPlinker_net.add_edge(type1_names[i],
                                  type2_names[i],
                                  weight=float(link_candidates[scorestr][i]))

        # Add edges between molecular family members
        if link_type == 'spec-gcf':
            type1_unique = (np.unique(link_candidates[type1str])).astype(int)
            map_spec_fam = data_links.mapping_spec["fam-id"]
            for type1 in type1_unique:
                if map_spec_fam[type1] > 0:  # if no singleton
                    members = data_links.family_members[int(
                        map_spec_fam[type1])][0]

                    # select other family members if among link candidates
                    members_present = [
                        x for x in list(members) if x in list(type1_unique)
                    ]
                    for member in members_present:
                        NPlinker_net.add_edge(type1str + str(type1),
                                              type1str + str(member))

        # export graph for drawing (e.g. using Cytoscape)
        nx.write_graphml(NPlinker_net, network_filename)

    def plot_candidates(self,
                        P_cutoff=0.8,
                        score_type='likescore',
                        score_cutoff=0,
                        type='mf-gcf'):
        """
        Plot best rated correlations between gcfs and spectra/families
        plot in form of seaborn clustermap
        """

        # Select score type
        if score_type == 'metcalf':
            scorestr = 'metcalf score'
        elif score_type == 'likescore':
            scorestr = 'likelihood score'
        else:
            raise Exception(
                "Wrong score_type given. Must be one of: 'metcalf', 'likescore' ."
            )

        if type == 'spec-gcf':
            link_candidates = self.link_candidates_gcf_spec
            selected_ids = np.where(
                (link_candidates["P(gcf|spec)"] > P_cutoff)
                & (link_candidates[scorestr] > score_cutoff))[0]
        elif type == 'mf-gcf':
            link_candidates = self.link_candidates_gcf_fam
            selected_ids = np.where(
                (link_candidates["P(gcf|fam)"] > P_cutoff)
                & (link_candidates[scorestr] > score_cutoff))[0]
        else:
            raise Exception("Wrong correlation 'type' given.")

        mapping_fams = np.unique(link_candidates.iloc[selected_ids, 0])
        mapping_gcfs = np.unique(link_candidates.iloc[selected_ids, 1])
        unique_fams = len(mapping_fams)
        unique_gcfs = len(mapping_gcfs)

        M_links = np.zeros((unique_fams, unique_gcfs))

        # define colors for different BGC classes
        bigscape_classes_dict = {
            "Others": "C0",
            "NRPS": "C1",
            "PKS-NRP_Hybrids": "C2",
            "PKSother": "C3",
            "PKSI": "C4",
            "RiPPs": "C5",
            "Saccharides": "C6",
            "Terpene": "C7"
        }
        col_colors = []

        # Create matrix with relevant link scores
        # TODO replace by better numpy method...
        for i in range(0, link_candidates.shape[0]):
            x = np.where(mapping_fams == link_candidates.iloc[i, 0])[0]
            y = np.where(mapping_gcfs == link_candidates.iloc[i, 1])[0]
            M_links[x, y] = link_candidates[scorestr][i]

        # make pandas dataframe from numpy array
        M_links = pd.DataFrame(M_links,
                               index=mapping_fams.astype(int),
                               columns=mapping_gcfs.astype(int))
        if type == 'spec-gcf':
            M_links.index.name = 'spectrum number'
        elif type == 'mf-gcf':
            M_links.index.name = 'molecular family number'
        M_links.columns.name = 'gene cluster family (GCF)'

        # add color label representing gene cluster class
        for bgc_class in link_candidates["BGC class"]:
            if bgc_class is None:
                col_colors.append((0, 0, 0))
            else:
                col_colors.append(bigscape_classes_dict[bgc_class])

    #    bgc_type_colors = pd.Series(mapping_gcfs.astype(int), index=M_links.columns).map(col_colors)
        graph = sns.clustermap(
            M_links,
            metric="correlation",
            method="weighted",
            cmap=
            "Reds",  #sns.cubehelix_palette(8, start=0, rot=.3, dark=0, light=1),
            vmin=0,
            vmax=np.max(link_candidates[scorestr]),
            col_cluster=True,
            col_colors=col_colors,
            robust=True)
        graph.fig.suptitle('Correlation map')

        # Rotate labels
        plt.setp(graph.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
        plt.setp(graph.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)

        # Make labels smaller
        plt.setp(graph.ax_heatmap.xaxis.get_majorticklabels(), fontsize=7)
        plt.setp(graph.ax_heatmap.yaxis.get_majorticklabels(), fontsize=7)

        plt.ylabel("scoring index")

        return M_links
