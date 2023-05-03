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

# Some functions for processing results - should be put somewhere else
import numpy as np


def get_sig_links(scores,
                  random_scores,
                  p_threshold=0.95,
                  direction='greater'):
    # gcfs are columns
    n_spec, n_gcf = scores.shape
    sig_links = np.zeros((n_spec, n_gcf))
    for gpos in range(n_gcf):
        max_rand = random_scores[:, gpos].max()
        min_rand = random_scores[:, gpos].min()
        if p_threshold > 1:
            for i, s in enumerate(scores[:, gpos]):
                if direction == 'greater':
                    if s > max_rand:
                        sig_links[i, gpos] = 1
                else:
                    if s < min_rand:
                        sig_links[i, gpos] = 1
        else:
            perc = np.percentile(random_scores[:, gpos],
                                 int(100 * p_threshold))
            for i, s in enumerate(scores[:, gpos]):
                if direction == 'greater':
                    if s >= perc:
                        sig_links[i, gpos] = 1
                else:
                    if s <= perc:
                        sig_links[i, gpos] = 1
    return sig_links


def get_sig_spec(data_link, sig_links, scores, gcf_pos, min_n_strains=2):
    # Check if there are *any* strains in the GCF
    # No strains = MiBIG
    # Can also filter if only (e.g. 2 strains)
    strain_sum = data_link.occurrence_gcf_strain[gcf_pos, :].sum()
    if strain_sum < min_n_strains:
        return []
    col = sig_links[:, gcf_pos]  # get the column
    sig_pos = np.where(col == 1)[0]
    orig_ids = []
    for sp in sig_pos:
        orig_ids.append(
            (int(data_link.mapping_spec.iloc[sp]["original spec-id"]),
             scores[sp, gcf_pos]))
    orig_ids.sort(key=lambda x: x[1], reverse=True)
    return orig_ids
