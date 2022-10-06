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
from nplinker.annotations import GNPS_DATA_COLUMNS
from nplinker.annotations import GNPS_KEY
from nplinker.annotations import create_gnps_annotation
from nplinker.logconfig import LogConfig
from nplinker.parsers.mgf import LoadMGF
from .load_gnps import load_gnps
from .MolecularFamily import MolecularFamily
from .SingletonFamily import SingletonFamily
from .Spectrum import Spectrum


logger = LogConfig.getLogger(__file__)


def mols_to_spectra(ms2, metadata):
    ms2_dict = {}
    for m in ms2:
        if not m[3] in ms2_dict:
            ms2_dict[m[3]] = []
        ms2_dict[m[3]].append((m[0], m[2]))

    spectra = []
    for i, m in enumerate(ms2_dict.keys()):
        new_spectrum = Spectrum(i, ms2_dict[m], int(m.name),
                                metadata[m.name]['precursormass'],
                                metadata[m.name]['parentmass'])
        new_spectrum.metadata = metadata[m.name]
        # add GNPS ID if in metadata under SPECTRUMID (this doesn't seem to be in regular MGF files
        # from GNPS, but *is* in the rosetta mibig MGF)
        # note: LoadMGF seems to lowercase (some) metadata keys?
        if 'spectrumid' in new_spectrum.metadata:
            # add an annotation for consistency, if not already there
            if GNPS_KEY not in new_spectrum.annotations:
                gnps_anno = {k: None for k in GNPS_DATA_COLUMNS}
                gnps_anno['SpectrumID'] = new_spectrum.metadata['spectrumid']
                create_gnps_annotation(new_spectrum, gnps_anno)
        spectra.append(new_spectrum)

    return spectra


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
    # common steps to all formats of GNPS data:
    #   - parse the MGF file to create a set of Spectrum objects
    #   - parse the edges file and update the spectra with that data

    # build a set of Spectrum objects by parsing the MGF file
    spec_dict = load_spectra(mgf_file, edges_file)

    unknown_strains = load_gnps(strains, nodes_file, quant_table_file,
                                metadata_table_file, ext_metadata_parsing,
                                spec_dict)

    molfams = _make_families(spec_dict.values())
    return spec_dict, list(spec_dict.values()), molfams, unknown_strains


def load_spectra(mgf_file, edges_file):
    ms1, ms2, metadata = LoadMGF(name_field='scans').load_spectra([mgf_file])
    logger.info(f'{len(ms1)} molecules parsed from MGF file')
    spectra = mols_to_spectra(ms2, metadata)
    # above returns a list, create a dict indexed by spectrum_id to make
    # the rest of the parsing a bit simpler
    spec_dict = {spec.spectrum_id: spec for spec in spectra}

    # add edges info to the spectra
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
