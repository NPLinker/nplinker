from __future__ import annotations
from deprecated import deprecated
from nplinker.logconfig import LogConfig
from nplinker.scoring.linking.utils import calc_likelihood_matrix


logger = LogConfig.getLogger(__file__)


@deprecated(version="1.3.3", reason="It's unused and will be removed in 2.0.0")
class LinkLikelihood:
    """Class to:
    1) create ansd store likelihood matrices (from co-occurences)
    2) select potential calculates for links.
    """

    def __init__(self):
        """Matrices that store likelihoods of empirically found co-occurence
        Example:
            P_spec_givengcf contains likelihoods P(spec_x|gcf_y),
            which is the probability of finding spec_x given there is a gcf_y.
        """
        # Probabilities reflecting co-occurences spectra <-> GCFs
        self.P_spec_given_gcf = []
        self.P_spec_not_gcf = []
        self.P_gcf_given_spec = []
        self.P_gcf_not_spec = []
        # Probabilities reflecting co-occurences mol.families <-> GCFs
        self.P_fam_given_gcf = []
        self.P_fam_not_gcf = []
        self.P_gcf_given_fam = []
        self.P_gcf_not_fam = []

    def calculate_likelihoods(self, data_links, type="spec-gcf"):
        """Calulate likelihoods from empirically found co-occurences in data.

        IF type='spec-gcf':
        P(GCF_x | spec_y), P(spec_y | GCF_x),
        P(GCF_x | not spec_y), P(spec_y | not GCF_x)
        IF type='mf-gcf':
        P(GCF_x | fam_y), P(fam_y | GCF_x),
        P(GCF_x | not fam_y), P(fam_y | not GCF_x)
        """
        # Make selection for scenario spec<->gcf or mf<->gcf
        if type == "spec-gcf":
            M_type1_type2 = data_links.cooccurrence_spec_gcf
            M_type1_nottype2 = data_links.cooccurrence_spec_notgcf
            M_nottype1_type2 = data_links.cooccurrence_notspec_gcf
            M_type1_cond = data_links.occurrence_spec_strain
        elif type == "mf-gcf":
            M_type1_type2 = data_links.cooccurrence_mf_gcf
            M_type1_nottype2 = data_links.cooccurrence_mf_notgcf
            M_nottype1_type2 = data_links.cooccurrence_notmf_gcf
            M_type1_cond = data_links.occurrence_mf_strain
        elif type == "spec-bgc" or type == "fam-bgc":
            raise Exception("Given types are not yet supported... ")
        else:
            raise Exception(
                "Wrong correlation 'type' given. Must be one of 'spec-gcf', 'mf-gcf'..."
            )

        logger.debug(f"Calculating likelihood matrices of type: {type}")
        # Calculate likelihood matrices using calc_likelihood_matrix()
        (
            P_type2_given_type1,
            P_type2_not_type1,
            P_type1_given_type2,
            P_type1_not_type2,
        ) = calc_likelihood_matrix(
            M_type1_cond,
            data_links.occurrence_gcf_strain,
            M_type1_type2,
            M_type1_nottype2,
            M_nottype1_type2,
        )
        if type == "spec-gcf":
            self.P_gcf_given_spec = P_type2_given_type1
            self.P_gcf_not_spec = P_type2_not_type1
            self.P_spec_given_gcf = P_type1_given_type2
            self.P_spec_not_gcf = P_type1_not_type2
        elif type == "mf-gcf":
            self.P_gcf_given_fam = P_type2_given_type1
            self.P_gcf_not_fam = P_type2_not_type1
            self.P_fam_given_gcf = P_type1_given_type2
            self.P_fam_not_gcf = P_type1_not_type2
        else:
            raise Exception("No correct likelihood matrices were created.")
