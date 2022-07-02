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

import os
from Bio import SeqIO
from .logconfig import LogConfig


logger = LogConfig.getLogger(__file__)


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
