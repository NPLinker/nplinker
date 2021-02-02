import os
import re

from ..logconfig import LogConfig
logger = LogConfig.getLogger(__file__)

# this will match strings like '...cluster001.gbk' or '...region022.gbk',
# and allow the number to be extracted easily
CLUSTER_REGION_REGEX = re.compile('(.+?)\\.(cluster|region)(\\d+).gbk$')

class KCBParser(object):
    """
    Parser for antismash knownclusterblast text output files
    """

    def __init__(self, filename):
        if not os.path.exists(filename):
            raise Exception('KCBParser failed to find file "{}"'.format(filename))

        self.bgc_genes = set()
        self.mibig_bgcs = []
        self.hits = {}

        with open(filename, 'r') as f:
            line = next(f)
            while not line.startswith('Table of genes'):
                line = next(f)
            # now we're in the first block. this is the section beginning with the
            # line "Table of genes, locations, strands and annotations of query cluster:"
            # and ending at "Significant hits:"
            top_block = []
            while True:
                line = next(f)
                if line.startswith('Significant'):
                    break
                else:
                    if len(line) > 1:
                        top_block.append(line.rstrip())

            # now we're in the second block. this is the section beginning with "Significant hits:"
            # and ending at "Details:"
            second_block = []
            while True:
                line = next(f)
                if line.startswith('Details'):
                    break
                else:
                    if len(line) > 1:
                        second_block.append(line.rstrip())


            # find the start of the first ">>" delimited section
            while True:
                try:
                    line = next(f)
                    if line.startswith('>>'):
                        break
                except Exception as e:
                    # EOF
                    return

            details = []
            finished = False
            while not finished:
                temp_list = []
                while True:
                    try:
                        line = next(f)
                        if line.startswith('>>'):
                            # finished current section
                            details.append(temp_list)
                            break
                        else:
                            if len(line) > 1:
                                temp_list.append(line.rstrip())
                    except Exception as e:
                        details.append(temp_list)
                        finished = True
                        break

            # do some processing on the blocks
            # firstly, extract the genes from the BGC -- stored in the first block
            # these are the gene names in the BGC
            for line in top_block:
                tokens = line.split()
                self.bgc_genes.add(tokens[0])

            # secondly, extract the MiBIG BGCs that are mentioned here
            for line in second_block:
                tokens = line.split()
                bgc_id = tokens[1]
                bgc_product_name = tokens[2]
                self.mibig_bgcs.append((bgc_id, bgc_product_name))

            for i, detail in enumerate(details):
                current_bgc_id = detail[0].split()[1] # this is the MiBIG ID, e.g. BGC0001666
                self.hits[current_bgc_id] = {}
                self.hits[current_bgc_id]['all_bgc_genes'] = self.bgc_genes

                assert current_bgc_id == self.mibig_bgcs[i][0] # they should be in the same order
                table_pos = detail.index('Table of genes, locations, strands and annotations of subject cluster:')
                pos = detail.index('Table of Blast hits (query gene, subject gene, %identity, blast score, %coverage, e-value):')

                # get the list of all of the genes within the MiBIG BGC
                all_mibig_genes = []
                for line in detail[table_pos+1:pos]:
                    all_mibig_genes.append(line.split()[0])

                self.hits[current_bgc_id]['all_mibig_genes'] = all_mibig_genes
                self.hits[current_bgc_id]['individual_hits'] = []

                for line in detail[pos+1:]:
                    tokens = line.split()
                    bgc_id = tokens[0]
                    self.hits[current_bgc_id]['individual_hits'].append({'source_bgc_gene': tokens[0], 
                                                      'mibig_bgc_gene': tokens[1],
                                                      'identity_percent': int(tokens[2]),
                                                      'blast_score': int(tokens[3])})

    @staticmethod
    def get_kcb_filename_from_bgc(bgc):
        """Given a BGC object, return the filename of the corresponding knownclusterblast .txt file (if any).

            
        This method attempts to derive the name of the knownclusterblast output file for the 
        supplied BGC object, using the original path + filename of the .gbk that the BGC was 
        sourced from during the loading process. 

        Assumptions made:
            - the .gbk files are grouped in subdirectories with <dataset>/antismash/
            - each of these subdirectories contains a "knownclusterblast" subdirectory
            - within that, there is a single .txt file for each .gbk 
        
        For some reason there doesn't appear to be a simple way of matching up the .gbk/.txt
        files using their contents, only the filenames. And there are various different 
        possibilities for the naming schemes, so this method tries to account for the
        known variants discovered so far. 

        Args:
            bgc: a BGC object

        Returns:
            A string containing the filename of the knownclusterblast .txt file for the
            supplied BGC, or None if an error occurred/file doesn't exist
        """
        if bgc.antismash_file is None:
            logger.warning('BGC {} has no antismash_file set'.format(bgc))
            return None

        # expecting to find the .txt files inside a 'knownclusterblast' subdir in the 
        # same location as the .gbk file, give up if that doesn't exist
        base_path = os.path.join(os.path.dirname(bgc.antismash_file), 'knownclusterblast')

        if not os.path.exists(base_path):
            logger.warning('Expected "knownclusterblast" directory not found at "{}"'.format(base_path))
            return None

        # get the name of the .gbk file itself (no path)
        genbank_file = os.path.split(bgc.antismash_file)[1]

        if 'region' in genbank_file:
            # case 1: assume the genbank files have a <someID>.region<num>.gbk naming scheme.
            # this seems to be the case for all the more recent datasets i've seen.

            # use a regex to extract the portions of the filename we're interested in
            regex_obj = CLUSTER_REGION_REGEX.search(genbank_file)
            number = int(regex_obj.group(3))
            prefix = regex_obj.group(1)

            # construct the expected filename using the prefix and number
            # (note no leading zeroes on the number)
            kcb_name = os.path.join(base_path, '{}_c{}.txt'.format(prefix, number)) 
        elif 'cluster' in genbank_file:
            # case 2: assume the genbank files have a <someID>.cluster<num>.gbk naming scheme.
            # this is the case with the Crusemann dataset on the paired platform among others.

            # use a regex to extract the number from the filename
            regex_obj = CLUSTER_REGION_REGEX.search(genbank_file)
            number = int(regex_obj.group(3))

            # construct the expected filename
            # (note no leading zeroes on the number)
            kcb_name = os.path.join(base_path, 'cluster{}.txt'.format(number))
        else:
            logger.warning('Unknown GenBank file naming scheme, failed to determine knownclusterblast filenames!')
            return None

        return kcb_name
