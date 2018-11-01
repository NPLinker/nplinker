import argparse
import os

from Bio import SeqIO

from aalist import AA_LIST as AA_CODES


AA_CODES_ISO = {}
for c in AA_CODES:
    AA_CODES_ISO['d-%s' % c] = c


class AntiSmashFile(object):
    def __init__(self, filename):
        self.raw_data = []
        self.filename = filename

        with open(filename, 'rU') as f:
            for record in SeqIO.parse(f, 'gb'):
                self.raw_data.append(AntiSmashRecord(record))

    def get_spec(self):
        for r in self.raw_data:
            for s in r.get_spec():
                yield s

    @property
    def products(self):
        return [x.products for x in self.raw_data]

    def build_prob(self):
        [x.build_prob() for x in self.raw_data]

    def get_prob(self, aa):
        return [x.get_prob(aa) for x in self.raw_data]


class AntiSmashRecord(object):
    def __init__(self, seq_record):
        self.raw_data = seq_record
        self.description = self.raw_data.description

        clusters = [x for x in self.raw_data.features if x.type == 'cluster']

        products = [x.qualifiers['product'] for x in clusters]
        smiles = [[y for y in x.qualifiers['note'] if y.startswith('Monomers prediction')] for x in clusters]

        asdomains = [x for x in self.raw_data.features if x.type == 'aSDomain']
        asdomains_with_predictions = [x for x in asdomains if 'specificity' in x.qualifiers]
        asdomains_with_predictions_known = [x for x in asdomains_with_predictions if 'AMP-binding' in x.qualifiers['domain']]

        self.products = products
        self.smiles = smiles
        self.asdomains_with_predictions_known = asdomains_with_predictions_known

        self.build_prob()

    def display(self):
        print(self.products)
        print(self.smiles)
        print(self.asdomains_with_predictions_known)

        print(self.asdomains_with_predictions_known[0].qualifiers['specificity'])

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
        self.specificities = [process_specificity(x.qualifiers['specificity']) for x in self.asdomains_with_predictions_known]

    def get_prob(self, aa):
        prob = 1.0
        for domain in self.specificities:
            domain_vote_sum = 0
            for vote in domain:
                if 'none' in vote:
                    p_vote = 1.0 / len(AA_CODES)
                elif aa in vote:
                    p_vote = 1.0 / len(vote)
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
        if prediction_string.startswith('PID to NN:') or prediction_string.startswith('SNN score:'):
            continue
        if prediction_string.startswith('Stachelhaus code:'):
            prediction = stachelhaus(prediction_string)
        elif prediction_string.startswith('NRPSpredictor3 SVM'):
            prediction = nrpspredictor3(prediction_string)
        elif prediction_string.startswith('pHMM'):
            prediction = phmm(prediction_string)
        elif prediction_string.startswith('PrediCAT'):
            prediction = predicat(prediction_string)
        elif prediction_string.startswith('SANDPUMA'):
            prediction = sandpuma(prediction_string)
        else:
            prediction = []
        votes.append(prediction)
    return votes


def to_set(nested):
    for item in nested:
        if type(item) is str:
            yield item
        else:
            for x in to_set(item):
                yield x


def predict_aa(filename):
    asr = AntiSmashFile(filename)
    for item in asr.raw_data:
        # for aa in item.get_probable_aa():
        for aa in AA_CODES:
            res = item.get_prob(aa)
            yield aa, res


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Predict AA specificity for gbk file')
    parser.add_argument('file', help='.gbk file with specificity predictions')
    args = parser.parse_args()

    for aa, res in predict_aa(args.file):
        print('%s,%s' % (aa, res))
