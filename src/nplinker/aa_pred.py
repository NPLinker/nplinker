# Copyright 2021 The NPLinker Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
from Bio import SeqIO
from .aalist import AA_LIST as AA_CODES


# AA_CODES_ISO contains codes that should map to
# an entry in the AA_CODES list - i.e. isomers, etc.
AA_CODES_ISO = {}
for c in AA_CODES:
    AA_CODES_ISO['d-%s' % c] = c
AA_CODES_ISO['b-ala'] = 'ala'


class AntiSmashFile():

    def __init__(self, filename):
        self.raw_data = []
        self.filename = filename

        with open(filename) as f:
            # antiSMASH 5 vs. 4 output is vastly different
            for record in SeqIO.parse(f, 'gb'):
                if 'structured_comment' in record.annotations:
                    as_version = record.annotations['structured_comment'][
                        'antiSMASH-Data']['Version']
                    as_major_version = as_version.split('.')[0]
                    if as_major_version == '5':
                        self.raw_data.append(AntiSmash5Record(record))
                    else:
                        raise ValueError(
                            f'Invalid antiSMASH version: {as_version}')
                else:
                    self.raw_data.append(AntiSmashRecord(record))

    def get_spec(self):
        for r in self.raw_data:
            yield from r.get_spec()

    @property
    def products(self):
        return [x.products for x in self.raw_data]

    def build_prob(self):
        [x.build_prob() for x in self.raw_data]

    def get_prob(self, aa):
        return [x.get_prob(aa) for x in self.raw_data]


class AntiSmash5Record():

    def __init__(self, seq_record):
        self.raw_data = seq_record
        self.description = self.raw_data.description

        clusters = [
            x for x in self.raw_data.features if x.type == 'cand_cluster'
        ]

        products = []
        for prod_list in [x.qualifiers['product'] for x in clusters]:
            for prod in prod_list:
                products.append(prod)

        smiles = []
        for smiles_list in [x.qualifiers['SMILES'] for x in clusters]:
            for s in smiles_list:
                smiles.append(s)

        asdomains = [x for x in self.raw_data.features if x.type == 'aSDomain']
        asdomains_with_predictions = [
            x for x in asdomains if 'specificity' in x.qualifiers
        ]
        asdomains_with_predictions_known = [
            x for x in asdomains_with_predictions
            if 'AMP-binding' in x.qualifiers['aSDomain']
        ]

        self.products = products
        self.smiles = smiles
        self.asdomains_with_predictions_known = asdomains_with_predictions_known

        self.specificities = []
        specificities_listed = [
            x.qualifiers['specificity']
            for x in asdomains_with_predictions_known
        ]
        for specificities_single_list in specificities_listed:
            for specificity in specificities_single_list:
                if specificity.startswith('consensus: '):
                    self.specificities.append(specificity.split()[1])

    def get_prob(self, aa):
        if aa in self.specificities:
            return 1.0
        else:
            return 0.0

    def get_spec(self):
        for aa in set(self.specificities):
            if aa in AA_CODES:
                yield aa
            elif aa in AA_CODES_ISO:
                yield AA_CODES_ISO[aa]


class AntiSmashRecord():

    def __init__(self, seq_record):
        self.raw_data = seq_record
        self.description = self.raw_data.description

        clusters = [x for x in self.raw_data.features if x.type == 'cluster']

        products = [x.qualifiers['product'] for x in clusters]
        smiles = [[
            y for y in x.qualifiers['note']
            if y.startswith('Monomers prediction')
        ] for x in clusters]

        asdomains = [x for x in self.raw_data.features if x.type == 'aSDomain']
        asdomains_with_predictions = [
            x for x in asdomains if 'specificity' in x.qualifiers
        ]
        asdomains_with_predictions_known = [
            x for x in asdomains_with_predictions
            if 'AMP-binding' in x.qualifiers['domain']
        ]

        self.products = products
        self.smiles = smiles
        self.asdomains_with_predictions_known = asdomains_with_predictions_known

        self.build_prob()

    def display(self):
        print(self.products)
        print(self.smiles)
        print(self.asdomains_with_predictions_known)

        print(
            self.asdomains_with_predictions_known[0].qualifiers['specificity'])

    def get_probable_aa(self):
        aalist = []
        for domain in self.specificities:
            for vote in domain:
                aalist.extend(vote)
        aaset = set(aalist)
        if 'none' in aaset:
            aaset.remove('none')
        return aaset

    def get_spec(self):
        for domain in self.asdomains_with_predictions_known:
            yield process_specificity(domain.qualifiers['specificity'])
            # for x in domain.qualifiers['specificity']:
            #     if x.startswith('PID to NN:') or x.startswith('SNN score:'):
            #         continue
            #     else:
            #         yield x

            #x_split = x.split(': ')
            #if len(x_split) > 1:
            #    yield x_split[1]
            #else:
            #    x_split = x.split(' ')
            #    yield x_split[1]

    def build_prob(self):
        self.specificities = [
            process_specificity(x.qualifiers['specificity'])
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
                    if vote == 'none':
                        p_vote = 0
                    #     p_vote = 1.0 / len(AA_CODES)
                    elif vote == aa:
                        p_vote = 1.0 / len(set(votes))
                    else:
                        p_vote = 0
                    domain_vote_sum += p_vote
            domain_prob = domain_vote_sum / len(domain)
            prob *= (1 - domain_prob)

        return 1 - prob


def stachelhaus(prediction_string):
    prediction_string = prediction_string.split(':')[-1].strip()
    predictions = predict_aa_string(prediction_string)
    return predictions


def sandpuma(prediction_string):
    prediction_string = prediction_string.split(':')[-1].strip()
    predictions = predict_aa_string(prediction_string)
    return predictions


def predicat(prediction_string):
    prediction_string = prediction_string.split('-')[-1]
    predictions = predict_aa_string(prediction_string)
    return predictions


def phmm(prediction_string):
    prediction_string = prediction_string.split(':')[-1].strip()
    predictions = predict_aa_string(prediction_string)
    return predictions


def nrpspredictor3(prediction_string):
    prediction_string = prediction_string.split(':')[-1].strip()
    predictions = predict_aa_string(prediction_string)
    return predictions


def predict_aa_string(prediction_string):
    parsed_predictions = []
    for prediction in prediction_string.split('|'):
        if prediction in AA_CODES:
            parsed_predictions.append(prediction)
        elif prediction in AA_CODES_ISO:
            parsed_predictions.append(AA_CODES_ISO[prediction])
        else:
            parsed_predictions.append('none')
    return parsed_predictions


def process_specificity(prediction_list):
    votes = []
    for prediction_string in prediction_list:
        if prediction_string.startswith(
                'PID to NN:') or prediction_string.startswith('SNN score:'):
            continue
        if prediction_string.startswith('SANDPUMA'):
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
    """
    Read AA losses from data file. (assume fixed structure...)
    """
    aa_losses = {}
    with open(filename) as f:
        reader = csv.reader(f, delimiter=',')
        next(reader)  # skip headers
        for line in reader:
            if len(line) == 0:
                continue
            aa_id = line[1]
            aa_mono = float(line[4])
            aa_avg = float(line[5])
            aa_losses[aa_id.lower()] = (aa_mono, aa_avg)

    return aa_losses


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Predict AA specificity for gbk file')
    parser.add_argument('file', help='.gbk file with specificity predictions')
    args = parser.parse_args()

    for aa, res in predict_aa(args.file):
        print(f'{aa},{res}')
