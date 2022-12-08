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

import csv
from os import PathLike

from nplinker.logconfig import LogConfig
from nplinker.metabolomics.gnps.gnps_spectrum_loader import GNPSSpectrumLoader
from nplinker.metabolomics.spectrum import Spectrum
from .load_gnps import load_gnps
from .molecular_family import MolecularFamily
from .singleton_family import SingletonFamily


logger = LogConfig.getLogger(__name__)

def load_edges(edges_file, spec_dict):
    logger.debug('loading edges file: {} [{} spectra from MGF]'.format(
        edges_file, len(spec_dict)))
    with open(edges_file) as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        try:
            cid1_index = headers.index('CLUSTERID1')
            cid2_index = headers.index('CLUSTERID2')
            cos_index = headers.index('Cosine')
            fam_index = headers.index('ComponentIndex')
        except ValueError as ve:
            raise Exception(
                'Unknown or missing column(s) in edges file: {}'.format(
                    edges_file))

        for line in reader:
            spec1_id = int(line[cid1_index])
            spec2_id = int(line[cid2_index])
            cosine = float(line[cos_index])
            family = int(line[fam_index])

            spec1 = spec_dict[spec1_id]
            spec2 = spec_dict[spec2_id]

            if family != -1:  # singletons
                spec1.family_id = family
                spec2.family_id = family

                spec1.edges.append((spec2.id, spec2.spectrum_id, cosine))
                spec2.edges.append((spec1.id, spec1.spectrum_id, cosine))
            else:
                spec1.family_id = family


def load_dataset(strains,
                 mgf_file,
                 edges_file,
                 nodes_file,
                 quant_table_file=None,
                 metadata_table_file=None,
                 ext_metadata_parsing=False):
    """Wrapper function instead of GNPSLoader.
    Loads all data required for the metabolomics side.

    Args:
        strains(StrainCollection): Object holding the strain collections
        mgf_file(string): Path to spectra file in MGF format
        edges_file(string): Path to GNPS edges.pairinfo file.
        nodes_file(string): Path to GNPS networking output file in tsv format.
        quant_table_file (string, optional): Path to quantification/annotation table. Defaults to None.
        metadata_table_file (string, optional): Unknown. Defaults to None.
        ext_metadata_parsing (bool, optional): Unknown. Defaults to False.

    Returns:
        Tuple[dict, List[Spectrum], List[MolecularFamily], dict]: Loaded Metabolomics data.
    """
    # common steps to all formats of GNPS data:
    #   - parse the MGF file to create a set of Spectrum objects
    #   - parse the edges file and update the spectra with that data

    # build a set of Spectrum objects by parsing the MGF file
    spec_dict = load_spectra(mgf_file, edges_file)

    # spectra = GNPSSpectrumLoader(mgf_file).spectra()
    # above returns a list, create a dict indexed by spectrum_id to make
    # the rest of the parsing a bit simpler
    # spec_dict = {spec.spectrum_id: spec for spec in spectra}

    # add edges info to the spectra
    molfams = _make_families(spec_dict.values())
    # molfams = GNPSMolecularFamilyLoader(edges_file).families()

    unknown_strains = load_gnps(strains, nodes_file, quant_table_file,
                                metadata_table_file, ext_metadata_parsing,
                                spec_dict)

    return spec_dict, list(spec_dict.values()), molfams, unknown_strains


def load_spectra(mgf_file: str | PathLike, edges_file: str | PathLike) -> dict[int, Spectrum]:
    """Wrapper function to load spectra and init the molecular family links.

    Args:
        mgf_file(str | PathLike): File storing the mass spectra in MGF format.
        edges_file(str | PathLike): File storing the molecular family information in .selfloop or .pairsinfo format.

    Returns:
        dict[int, Spectrum]: Indexed dict of mass spectra.
    """

    spectra = GNPSSpectrumLoader(mgf_file).spectra()
    # above returns a list, create a dict indexed by spectrum_id to make
    # the rest of the parsing a bit simpler
    spec_dict: dict[int, Spectrum] = {spec.spectrum_id: spec for spec in spectra}
    load_edges(edges_file, spec_dict)
    return spec_dict


def _make_families(spectra):
    families = []
    family_dict = {}
    family_index = 0
    fams, singles = 0, 0
    for spectrum in spectra:
        family_id = spectrum.family_id
        if family_id == -1:  # singleton
            new_family = SingletonFamily()
            new_family.id = family_index
            family_index += 1

            new_family.add_spectrum(spectrum)
            spectrum.family = new_family
            families.append(new_family)
            singles += 1
        else:
            if family_id not in family_dict:
                new_family = MolecularFamily(family_id)
                new_family.id = family_index
                new_family.family_id = family_id  # preserve original ID here
                family_index += 1

                new_family.add_spectrum(spectrum)
                spectrum.family = new_family
                families.append(new_family)
                family_dict[family_id] = new_family
                fams += 1
            else:
                family_dict[family_id].add_spectrum(spectrum)
                spectrum.family = family_dict[family_id]

    logger.debug('make_families: {} molams + {} singletons'.format(
        fams, singles))
    return families
