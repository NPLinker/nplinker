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
# CG: these functions are not used by other, could be removed
import logging
from . import nplinker_iokr
from .spectrum import MSSpectrum


logger = logging.getLogger(__name__)


def run_iokr_ranking(spec, bgc_list):
    logger.debug(f"IOKR datapath: {nplinker_iokr.get_datapath()}")

    # extract SMILES strings for each supplied BGC (caching them in the object)
    smiles = [bgc.smiles for bgc in bgc_list]

    # TODO could probably handle this more gracefully by not submitting any BGC with no SMILES
    if None in smiles:
        logger.error("Failed to run IOKR scoring, one or more BGCs missing SMILES data")
        return None

    iokr_server = nplinker_iokr.get_iokr_server()
    rank = iokr_server.rank_smiles(MSSpectrum(spec=spec), smiles)
    return [bgc_list[i] for i in rank]


def run_iokr_scoring(spec_list, bgc_list):
    logger.debug(f"IOKR datapath: {nplinker_iokr.get_datapath()}")
    smiles = [bgc.smiles for bgc in bgc_list]
    if None in smiles:
        logger.error("Failed to run IOKR scoring, one or more BGCs missing SMILES data")
        return None

    spectra = [MSSpectrum(spec=spec) for spec in spec_list]

    iokr_server = nplinker_iokr.get_iokr_server()
    rank = iokr_server.score_smiles(spectra, smiles)

    return rank
