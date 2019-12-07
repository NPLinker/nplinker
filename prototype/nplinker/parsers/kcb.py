import os

class KCBParser(object):
    """
    Parser for antismash knownclusterblast text output files
    """

    def __init__(self, filename):
        if not os.path.exists(filename):
            raise Exception('KCBParser: file not found: "{}"'.format(filename))

        self.bgc_genes = set()
        self.mibig_bgcs = []
        self.hits = {}

        with open(filename, 'r') as f:
            line = next(f)
            while not line.startswith('Table of genes'):
                line = next(f)
            # now we're in the first block
            top_block = []
            while True:
                line = next(f)
                if line.startswith('Significant'):
                    break
                else:
                    if len(line) > 1:
                        top_block.append(line.rstrip())
            # now we're in the second block
            second_block = []
            while True:
                line = next(f)
                if line.startswith('Details'):
                    break
                else:
                    if len(line) > 1:
                        second_block.append(line.rstrip())
            while True:
                try:
                    line = next(f)
                    if line.startswith('>>'):
                        break
                except:
                    return
            details = []
            finished = False
            while not finished:
                temp_list = []
                while True:
                    try:
                        line = next(f)
                        if line.startswith('>>'):
                            # finished one
                            details.append(temp_list)
                            break
                        else:
                            if len(line) > 1:
                                temp_list.append(line.rstrip())
                    except:
                        details.append(temp_list)
                        finished = True
                        break
            # do some processing on the blocks
            # firstly, extract the genes from the BGC -- stored in the first block
            for line in top_block:
                tokens = line.split()
                self.bgc_genes.add(tokens[0])
            # secondly, extract the BGCs that are mentioned here
            for line in second_block:
                tokens = line.split()
                bgc_id = tokens[1]
                bgc_product_name = tokens[2]
                self.mibig_bgcs.append((bgc_id,bgc_product_name))
            for i,detail in enumerate(details):
                current_bgc_id = detail[0].split()[1]
                self.hits[current_bgc_id] = []
                assert current_bgc_id == self.mibig_bgcs[i][0] # they should be in the same order
                table_pos = detail.index('Table of genes, locations, strands and annotations of subject cluster:')
                pos = detail.index('Table of Blast hits (query gene, subject gene, %identity, blast score, %coverage, e-value):')
                for line in detail[pos+1:]:
                    tokens = line.split()
                    bgc_id = tokens[0]

                    self.hits[current_bgc_id].append({'source_bgc_gene': tokens[0], 
                                                      'mibig_bgc_gene': tokens[1],
                                                      'identity_percent': int(tokens[2]),
                                                      'blast_score': int(tokens[3]),
                                                      'all_bgc_genes': self.bgc_genes})

    @staticmethod
    def get_kcb_filename_from_bgc(bgc):
        if bgc.antismash_file is None:
            return None

        base_path = os.path.join(os.path.dirname(bgc.antismash_file), 'knownclusterblast')
        genbank_file = bgc.antismash_file.split(os.sep)[-1]
        # remove the "regionXXX" and turn it into "_cXXX"
        tokens = genbank_file.split('region')
        number = int(tokens[1].split('.')[0])
        start_name = tokens[0][:-1] # remove the last dot
        start_name += '_c{}.txt'.format(number)
        kcb_name = os.path.join(base_path, start_name)
        return kcb_name

if __name__ == "__main__":
    kcbp = KCBParser('KRD012.Scaffold_20_c1.txt')
    print(kcbp.hits)

