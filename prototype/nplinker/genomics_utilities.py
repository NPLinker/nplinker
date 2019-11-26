import os

from Bio import SeqIO

from logconfig import LogConfig
logger = LogConfig.getLogger(__file__)

def get_known_cluster_blast(bgc):
    if hasattr(bgc,'antismash_file'):           
        bgc_file = bgc.antismash_file
        if not bgc_file:
            print("No antismash file in {}".format(bgc))
            return [] # Stops it being repeatedly called
        records = list(SeqIO.parse(bgc.antismash_file,'genbank'))
        hits = []
        for record in records:
            for feature in record.features:
                if feature.type == 'cluster':
                    bgc_list = feature.qualifiers.get('knownclusterblast',None)
                    if bgc_list:
                        for bgc_item in bgc_list:
                            tokens = bgc_item.split('\t')
                            bgc_name = tokens[0].split('.')[1].strip()
                            similarity = tokens[1].split('(')[1].split('%')[0]
                            hits.append((bgc_name,similarity))
        return hits

def get_smiles(bgc):
    if bgc.antismash_file is None:
        return None

    if not os.path.exists(bgc.antismash_file):
        logger.warn('Missing antismash_file: {}'.format(bgc.antismash_file))
        return None

    with open(bgc.antismash_file, 'r') as f:
        for rec in SeqIO.parse(f, 'gb'):
            # rec is a Bio.SeqRecord object, search its .features list 
            # for "cand_cluster" and then extract SMILES string from there
            # TODO is this always correct or can it appear in other places?
            for feature in rec.features:
                if feature.type == 'cand_cluster':
                    smiles = feature.qualifiers.get('SMILES', None)
                    if smiles is None: 
                        return None

                    # always a list (TODO?)
                    if len(smiles[0]) == 0:
                        return None

                    # seem to get space chars in some of these, which are not allowed
                    # by the SMILES spec, so strip them out here
                    return smiles[0].replace(' ', '')
