# CG: This module parse AntiSMASH files to get the predicted monomers of product.
# And it only cares the case when the monomers are amino acids.
# TODO This module is used in bgc.py, referenced in gcf.py and then in misc.py,
# NPLinker core business does not use it at all.

from __future__ import annotations
import argparse
from Bio import SeqIO


AA_CODES = [
    "ala",
    "cys",
    "asp",
    "glu",
    "phe",
    "gly",
    "his",
    "ile",
    "lys",
    "leu",
    "met",
    "asn",
    "pyl",
    "pro",
    "gln",
    "arg",
    "ser",
    "thr",
    "sec",
    "val",
    "trp",
    "tyr",
]

# AA_CODES_ISOMER contains codes that should map to
# an entry in the AA_CODES list - i.e. isomers, etc.
AA_CODES_ISOMER = {}
for c in AA_CODES:
    AA_CODES_ISOMER["d-%s" % c] = c
AA_CODES_ISOMER["b-ala"] = "ala"


# CG: create one class for handling antismash output files and records
## since antismash output formats change along with version
## we should delegate all processes in other classes/functions to the methods of this class
## so that one type of data, one class covering all methodes needed
class AntiSmashFile:
    def __init__(self, filename):
        """Parse product and specificity (predicted monomer) from AntiSMASH file.

        This class is a wrapper for AntiSmash4Record and AntiSmash5Record.

        Deprecated:
            AntiSMASH file always contains one sequence record, so users should
            directly use AntiSmash5Record or AntiSmash4Record.

        Args:
            filename: AntiSMASH file path
        """
        self.raw_data = []
        self.filename = filename

        with open(filename) as f:
            # antiSMASH 5 vs. 4 output is vastly different
            # CG: TODO should use SeqIO.read since antismash file has only one sequence record
            for record in SeqIO.parse(f, "gb"):
                if "structured_comment" in record.annotations:
                    as_version = record.annotations["structured_comment"]["antiSMASH-Data"][
                        "Version"
                    ]
                    as_major_version = as_version.split(".")[0]
                    if as_major_version == "5":
                        self.raw_data.append(AntiSmash5Record(record))
                    else:
                        raise ValueError(f"Invalid antiSMASH version: {as_version}")
                else:
                    self.raw_data.append(AntiSmash4Record(record))

    def get_spec(self):
        """Get specificity (predicted monomer)."""
        for r in self.raw_data:
            yield from r.get_spec()

    @property
    def products(self):
        """Get list of products."""
        return [x.products for x in self.raw_data]

    def build_prob(self):
        [x.build_prob() for x in self.raw_data]

    def get_prob(self, aa):
        return [x.get_prob(aa) for x in self.raw_data]


class AntiSmash5Record:
    # CG: TODO input should also include antismash file
    def __init__(self, seq_record):
        """Parse product and specificity (predicted monomer of product) from
            AntiSMASH v5 data.

        Args:
            seq_record: SeqRecord of AntiSMASH
        """
        self.raw_data = seq_record
        self.description = self.raw_data.description

        # parse "cand_cluster"
        clusters = [x for x in self.raw_data.features if x.type == "cand_cluster"]

        products = []
        for prod_list in [x.qualifiers["product"] for x in clusters]:
            for prod in prod_list:
                products.append(prod)

        smiles = []
        for smiles_list in [x.qualifiers["SMILES"] for x in clusters if "SMILES" in x.qualifiers]:
            for s in smiles_list:
                smiles.append(s)

        # parse "aSDomain"
        asdomains = [x for x in self.raw_data.features if x.type == "aSDomain"]
        asdomains_with_predictions = [x for x in asdomains if "specificity" in x.qualifiers]
        asdomains_with_predictions_known = [
            x for x in asdomains_with_predictions if "AMP-binding" in x.qualifiers["aSDomain"]
        ]

        self.products = products
        self.smiles = smiles
        self.asdomains_with_predictions_known = asdomains_with_predictions_known

        # parse predicted monomer of product, i.e. "specificity"
        self.specificities = []
        specificities_listed = [
            x.qualifiers["specificity"] for x in asdomains_with_predictions_known
        ]
        for specificities_single_list in specificities_listed:
            for specificity in specificities_single_list:
                if specificity.startswith("consensus: "):
                    self.specificities.append(specificity.split()[1])

    def get_prob(self, aa):
        """Get probability of predicted amino acid.

        Args:
            aa: amino acid

        Returns:
            prediction probability
        """
        if aa in self.specificities:
            return 1.0
        else:
            return 0.0

    def get_spec(self):
        """Get predicted specificities (animo acids)."""
        for aa in set(self.specificities):
            if aa in AA_CODES:
                yield aa
            elif aa in AA_CODES_ISOMER:
                yield AA_CODES_ISOMER[aa]


class AntiSmash4Record:
    # CG: this is for antiSMASH version 4.0 results, it can be removed.

    def __init__(self, seq_record):
        self.raw_data = seq_record
        self.description = self.raw_data.description

        clusters = [x for x in self.raw_data.features if x.type == "cluster"]

        products = [x.qualifiers["product"] for x in clusters]
        smiles = [
            [y for y in x.qualifiers["note"] if y.startswith("Monomers prediction")]
            for x in clusters
        ]

        asdomains = [x for x in self.raw_data.features if x.type == "aSDomain"]
        asdomains_with_predictions = [x for x in asdomains if "specificity" in x.qualifiers]
        asdomains_with_predictions_known = [
            x for x in asdomains_with_predictions if "AMP-binding" in x.qualifiers["domain"]
        ]

        self.products = products
        self.smiles = smiles
        self.asdomains_with_predictions_known = asdomains_with_predictions_known

        self.build_prob()

    def display(self):
        print(self.products)
        print(self.smiles)
        print(self.asdomains_with_predictions_known)

        print(self.asdomains_with_predictions_known[0].qualifiers["specificity"])

    def get_probable_aa(self):
        aalist = []
        for domain in self.specificities:
            for vote in domain:
                aalist.extend(vote)
        aaset = set(aalist)
        if "none" in aaset:
            aaset.remove("none")
        return aaset

    def get_spec(self):
        for domain in self.asdomains_with_predictions_known:
            yield process_specificity(domain.qualifiers["specificity"])
            # for x in domain.qualifiers['specificity']:
            #     if x.startswith('PID to NN:') or x.startswith('SNN score:'):
            #         continue
            #     else:
            #         yield x

            # x_split = x.split(': ')
            # if len(x_split) > 1:
            #    yield x_split[1]
            # else:
            #    x_split = x.split(' ')
            #    yield x_split[1]

    def build_prob(self):
        self.specificities = [
            process_specificity(x.qualifiers["specificity"])
            for x in self.asdomains_with_predictions_known
        ]

    def get_prob(self, aa):
        # 'none' votes do not affect predictions, because then
        # all AAs are equally likely, so we can't use them to exclude.
        prob = 1.0
        for domain in self.specificities:
            domain_vote_sum = 0
            for votes in domain:
                for vote in set(votes):
                    if vote == "none":
                        p_vote = 0
                    #     p_vote = 1.0 / len(AA_CODES)
                    elif vote == aa:
                        p_vote = 1.0 / len(set(votes))
                    else:
                        p_vote = 0
                    domain_vote_sum += p_vote
            domain_prob = domain_vote_sum / len(domain)
            prob *= 1 - domain_prob

        return 1 - prob


def stachelhaus(prediction_string):
    prediction_string = prediction_string.split(":")[-1].strip()
    predictions = predict_aa_string(prediction_string)
    return predictions


def sandpuma(prediction_string):
    prediction_string = prediction_string.split(":")[-1].strip()
    predictions = predict_aa_string(prediction_string)
    return predictions


def predicat(prediction_string):
    prediction_string = prediction_string.split("-")[-1]
    predictions = predict_aa_string(prediction_string)
    return predictions


def phmm(prediction_string):
    prediction_string = prediction_string.split(":")[-1].strip()
    predictions = predict_aa_string(prediction_string)
    return predictions


def nrpspredictor3(prediction_string):
    prediction_string = prediction_string.split(":")[-1].strip()
    predictions = predict_aa_string(prediction_string)
    return predictions


def predict_aa_string(prediction_string):
    parsed_predictions = []
    for prediction in prediction_string.split("|"):
        if prediction in AA_CODES:
            parsed_predictions.append(prediction)
        elif prediction in AA_CODES_ISOMER:
            parsed_predictions.append(AA_CODES_ISOMER[prediction])
        else:
            parsed_predictions.append("none")
    return parsed_predictions


def process_specificity(prediction_list):
    votes = []
    for prediction_string in prediction_list:
        if prediction_string.startswith("PID to NN:") or prediction_string.startswith("SNN score:"):
            continue
        if prediction_string.startswith("SANDPUMA"):
            prediction = sandpuma(prediction_string)
        else:
            continue
        # Only process SANDPUMA votes for now (to not have to collate the votes)
        # if prediction_string.startswith('Stachelhaus code:'):
        #     prediction = stachelhaus(prediction_string)
        # elif prediction_string.startswith('NRPSpredictor3 SVM'):
        #     prediction = nrpspredictor3(prediction_string)
        # elif prediction_string.startswith('pHMM'):
        #     prediction = phmm(prediction_string)
        # elif prediction_string.startswith('PrediCAT'):
        #     prediction = predicat(prediction_string)
        # elif prediction_string.startswith('SANDPUMA'):
        #     prediction = sandpuma(prediction_string)
        # else:
        #     prediction = []
        votes.append(prediction)
    return votes


def to_set(nested):
    for item in nested:
        if type(item) is str:
            yield item
        else:
            yield from to_set(item)


def predict_aa(filename):
    asr = AntiSmashFile(filename)
    for item in asr.raw_data:
        # for aa in item.get_probable_aa():
        for aa in AA_CODES:
            res = item.get_prob(aa)
            yield aa, res


def read_aa_losses(filename):
    """Read AA losses from data file. (assume fixed structure...)."""
    aa_losses = {}
    with open(filename) as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)  # skip headers
        for line in reader:
            if len(line) == 0:
                continue
            aa_id = line[1]
            aa_mono = float(line[4])
            aa_avg = float(line[5])
            aa_losses[aa_id.lower()] = (aa_mono, aa_avg)

    return aa_losses


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Predict AA specificity for gbk file")
    parser.add_argument("file", help=".gbk file with specificity predictions")
    args = parser.parse_args()

    for aa, res in predict_aa(args.file):
        print(f"{aa},{res}")
