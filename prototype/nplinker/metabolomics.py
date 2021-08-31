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
import csv
import re

from .parsers.mgf import LoadMGF
from .strains import StrainCollection
from .utils import sqrt_normalise
from .annotations import GNPS_KEY, GNPS_DATA_COLUMNS, create_gnps_annotation

from .logconfig import LogConfig
logger = LogConfig.getLogger(__file__)

JCAMP = '##TITLE={}\\n' +\
        '##JCAMP-DX=nplinker vTODO\\n' +\
        '##DATA TYPE=Spectrum\\n' +\
        '##DATA CLASS=PEAKTABLE\\n' +\
        '##ORIGIN=TODO_DATASET_ID\\n' +\
        '##OWNER=nobody\\n' +\
        '##XUNITS=M/Z\\n' +\
        '##YUNITS=RELATIVE ABUNDANCE\\n' +\
        '##NPOINTS={}\\n' +\
        '##PEAK TABLE=(XY..XY)\\n' +\
        '{}\\n' +\
        '##END=\\n'

# compile a regex for matching .mzXML and .mzML strings 
RE_MZML_MZXML = re.compile('.mzXML|.mzML')

class Spectrum(object):
    
    def __init__(self, id, peaks, spectrum_id, precursor_mz, parent_mz=None, rt=None):
        self.id = id
        self.peaks = sorted(peaks, key=lambda x: x[0]) # ensure sorted by mz
        self.normalised_peaks = sqrt_normalise(self.peaks) # useful later
        self.n_peaks = len(self.peaks)
        self.max_ms2_intensity = max([intensity for mz, intensity in self.peaks])
        self.total_ms2_intensity = sum([intensity for mz, intensity in self.peaks])
        assert(isinstance(spectrum_id, int))
        self.spectrum_id = spectrum_id
        self.rt = rt
        self.precursor_mz = precursor_mz
        self.parent_mz = parent_mz
        self.gnps_id = None # CCMSLIB...
        # TODO should add intensity here too
        self.metadata = {}
        self.edges = []
        self.strains = StrainCollection()
        # this is a dict indexed by Strain objects (the strains found in this Spectrum), with 
        # the values being dicts of the form {growth_medium: peak intensity} for the parent strain
        self.growth_media = {}
        self.family_id = -1
        self.family = None
        # a dict indexed by filename, or "gnps" 
        self.annotations = {}
        self._losses = None
        self._jcamp = None

    def add_strain(self, strain, alias_or_name, growth_medium, peak_intensity):
        # adds the strain to the StrainCollection if not already there
        self.strains.add(strain)

        if alias_or_name not in self.growth_media:
            self.growth_media[alias_or_name] = {}

        if growth_medium is None:
            self.growth_media[alias_or_name].update({'unknown_medium_{}'.format(len(self.growth_media[alias_or_name])): peak_intensity})
            return

        if alias_or_name in self.growth_media and growth_medium in self.growth_media[alias_or_name]:
            raise Exception('Growth medium clash: {} / {} {}'.format(self, strain, growth_medium))
        
        self.growth_media[alias_or_name].update({growth_medium: peak_intensity})
        
    @property
    def is_library(self):
        return GNPS_KEY in self.annotations

    def set_annotations(self, key, data):
        self.annotations[key] = data

    @property
    def gnps_annotations(self):
        if GNPS_KEY not in self.annotations:
            return None

        return self.annotations[GNPS_KEY][0]

    def has_annotations(self):
        return len(self.annotations) > 0

    def get_metadata_value(self, key):
        val = self.metadata.get(key, None)
        return val

    def has_strain(self, strain):
        return strain in self.strains

    def get_growth_medium(self, strain):
        if strain not in self.strains:
            return None

        gms = self.growth_media[strain]
        return list(gms.keys())[0]

    def to_jcamp_str(self, force_refresh=False):
        if self._jcamp is not None and not force_refresh:
            return self._jcamp

        peakdata = '\\n'.join('{}, {}'.format(*p) for p in self.peaks)
        self._jcamp = JCAMP.format(str(self), self.n_peaks, peakdata)
        return self._jcamp

    def __str__(self):
        return "Spectrum(id={}, spectrum_id={}, strains={})".format(self.id, self.spectrum_id, len(self.strains))

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return self.id == other.id

    def __hash__(self):
        return self.id

    def __cmp__(self, other):
        if self.parent_mz >= other.parent_mz:
            return 1
        else:
            return -1

    def __lt__(self,other):
        if self.parent_mz <= other.parent_mz:
            return 1
        else:
            return 0

    # from molnet repo
    def keep_top_k(self, k=6, mz_range=50):
        # only keep peaks that are in the top k in += mz_range
        start_pos = 0
        new_peaks = []
        for mz, intensity in self.peaks:
            while self.peaks[start_pos][0] < mz - mz_range:
                start_pos += 1
            end_pos = start_pos

            n_bigger = 0
            while end_pos < len(self.peaks) and self.peaks[end_pos][0] <= mz + mz_range:
                if self.peaks[end_pos][1] > intensity:
                    n_bigger += 1
                end_pos += 1

            if n_bigger < k:
                new_peaks.append((mz, intensity))

        self.peaks = new_peaks
        self.n_peaks = len(self.peaks)
        if self.n_peaks > 0:
            self.normalised_peaks = sqrt_normalise(self.peaks)
            self.max_ms2_intensity = max([intensity for mz, intensity in self.peaks])
            self.total_ms2_intensity = sum([intensity for mz, intensity in self.peaks])
        else:
            self.normalised_peaks = []
            self.max_ms2_intensity = 0.0
            self.total_ms2_intensity = 0.0

    @property
    def losses(self):
        """
        All mass shifts in the spectrum, and the indices of the peaks
        """
        if self._losses is None:
            # populate loss table
            losses = []
            for i in range(len(self.peaks)):
                loss = self.precursor_mz - self.peaks[i][0]
                losses.append((loss, self.id, i))
# THIS SEEMED TO ME LIKE IT WOULD TAKE THE WRONG DIFFERENCES AS LOSSES:
# TODO: please check!
#                for j in range(i):
#                    loss = self.peaks[i][0] - self.peaks[j][0]
#                    losses.append((loss, i, j))
                
            # Sort by loss
            losses.sort(key=lambda x: x[0])
            self._losses = losses
        return self._losses

    def has_loss(self, mass, tol):
        """
        Check if the scan has the specified loss (within tolerance)
        """
        matched_losses = []

        idx = 0
        # Check losses in range [0, mass]
        while idx < len(self.losses) and self.losses[idx][0] <= mass:
            if mass - self.losses[idx][0] < tol:
                matched_losses.append(self.losses[idx])
            idx += 1

        # Add all losses in range [mass, mass+tol(
        while idx < len(self.losses) and self.losses[idx][0] < mass + tol:
            matched_losses.append(self.losses[idx])
            idx += 1

        return matched_losses

class MolecularFamily(object):
    def __init__(self, family_id):
        self.id = -1 
        self.family_id = family_id
        self.spectra = []
        self.family = None

    def has_strain(self, strain):
        for spectrum in self.spectra:
            if spectrum.has_strain(strain):
                return True

        return False

    @property
    def strains(self):
        strains = set([])
        for spectrum in self.spectra:
            strains = strains.union(spectrum.strains)
        return strains

    def add_spectrum(self, spectrum):
        self.spectra.append(spectrum)

    def __str__(self):
        return 'MolFam(family_id={}, spectra={})'.format(self.family_id, len(self.spectra))

    def __eq__(self, other):
        return self.id == other.id

    def __hash__(self):
        return self.id

class SingletonFamily(MolecularFamily):
    def __init__(self):
        super(SingletonFamily, self).__init__(-1)

    def __str__(self):
        return "Singleton molecular family"

#
# methods for parsing metabolomics data files
#

GNPS_FORMAT_UNKNOWN         = 'unknown'
GNPS_FORMAT_OLD_ALLFILES    = 'allfiles'
GNPS_FORMAT_OLD_UNIQUEFILES = 'uniquefiles'
GNPS_FORMAT_NEW_FBMN        = 'fbmn'

def get_headers(filename, delimiters=['\t', ',']):
    headers = None
    with open(filename, 'r') as f:
        headers = f.readline()
        for dl in delimiters:
            if len(headers.split(dl)) < 2:
                continue
            headers = headers.split(dl)
            break
    return headers

def identify_gnps_format(filename, has_quant_table):
    headers = get_headers(filename)
    
    if headers is None:
        return GNPS_FORMAT_UNKNOWN
    
    # first, check for AllFiles
    if 'AllFiles' in headers:
        # this should be an old-style dataset like Crusemann, with a single .tsv file
        # containing all the necessary info. The AllFiles column should contain pairs
        # of mzXML filenames and scan numbers in this format:
        #   filename1:scannumber1###filename2:scannumber2###... 
        return GNPS_FORMAT_OLD_ALLFILES
    elif 'UniqueFileSources' in headers:
        # this is a slightly newer-style dataset, e.g. MSV000084771 on the platform
        # it still just has a single .tsv file, but AllFiles is apparently replaced
        # with a UniqueFileSources column. There is also a UniqueFileSourcesCount 
        # column which just indicates the number of entries in the UniqueFileSources
        # column. If there are multiple entries the delimiter is a | character
        return GNPS_FORMAT_OLD_UNIQUEFILES
    elif has_quant_table:
        # if there is no AllFiles/UniqueFileSources, but we DO have a quantification
        # table file, that should indicate a new-style dataset like Carnegie
        # TODO check for the required header columns here too
        return GNPS_FORMAT_NEW_FBMN
    else:
        # if we don't match any of the above cases then it's not a recognised format
        return GNPS_FORMAT_UNKNOWN

def mols_to_spectra(ms2, metadata):
    ms2_dict = {}
    for m in ms2:
        if not m[3] in ms2_dict:
            ms2_dict[m[3]] = []
        ms2_dict[m[3]].append((m[0], m[2]))
    
    spectra = []
    for i, m in enumerate(ms2_dict.keys()):
        new_spectrum = Spectrum(i, ms2_dict[m], int(m.name), metadata[m.name]['precursormass'], metadata[m.name]['parentmass'])
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
    logger.debug('loading edges file: {} [{} spectra from MGF]'.format(edges_file, len(spec_dict)))
    with open(edges_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader) 
        try:
            cid1_index = headers.index('CLUSTERID1')
            cid2_index = headers.index('CLUSTERID2')
            cos_index = headers.index('Cosine')
            fam_index = headers.index('ComponentIndex')
        except ValueError as ve:
            raise Exception('Unknown or missing column(s) in edges file: {}'.format(edges_file))

        for line in reader:
            spec1_id = int(line[cid1_index])
            spec2_id = int(line[cid2_index])
            cosine = float(line[cos_index])
            family = int(line[fam_index])

            spec1 = spec_dict[spec1_id]
            spec2 = spec_dict[spec2_id]

            if family != -1: # singletons
                spec1.family_id = family
                spec2.family_id = family

                spec1.edges.append((spec2.id, spec2.spectrum_id, cosine))
                spec2.edges.append((spec1.id, spec1.spectrum_id, cosine))
            else:
                spec1.family_id = family

def messy_strain_naming_lookup(mzxml, strains):
    if mzxml in strains:
        # life is simple!
        return strains.lookup(mzxml)

    # 1. knock off the .mzXML and try again (using splitext should also handle
    # other extensions like .mzML here)
    mzxml_no_ext = os.path.splitext(mzxml)[0]
    if mzxml_no_ext in strains:
        return strains.lookup(mzxml_no_ext)

    # 2. if that doesn't work, try using everything up to the first -/_
    underscore_index = mzxml_no_ext.find('_')
    hyphen_index = mzxml_no_ext.find('-')
    mzxml_trunc_underscore = mzxml_no_ext if underscore_index == -1 else mzxml_no_ext[:underscore_index]
    mzxml_trunc_hyphen = mzxml_no_ext if hyphen_index == -1 else mzxml_no_ext[:hyphen_index]
    if underscore_index != -1 and mzxml_trunc_underscore in strains:
        return strains.lookup(mzxml_trunc_underscore)
    if hyphen_index != -1 and mzxml_trunc_hyphen in strains:
        return strains.lookup(mzxml_trunc_hyphen)

    # 3. in the case of original Crusemann dataset, many of the filenames seem to
    # match up to real strains with the initial "CN" missing ???
    for mzxml_trunc in [mzxml_trunc_hyphen, mzxml_trunc_underscore]:
        if 'CN' + mzxml_trunc in strains:
            return strains.lookup('CN' + mzxml_trunc)

    # give up
    return None

def md_convert(val):
    """Try to convert raw metadata values from text to integer, then float if that fails"""
    try:
        return int(val)
    except (ValueError, TypeError):
        try:
            return float(val)
        except (ValueError, TypeError):
            if val.lower() == 'n/a':
                return None
    return val

def parse_mzxml_header(hdr, strains, md_table, ext_metadata_parsing):
    # given a column header from the quantification_table file produced by GNPS,
    # this method checks if it's one that contains peak information for a given 
    # strain label by searching for the 'Peak area' string. 
    #    e.g. 'KRD168_ISP3.mzML Peak area'
    # it then extracts the strain label by taking the text up to the '.mzML' part
    # and attempts to match this to the set of strains provided by the user in
    # the strain_mappings file. finally it also tries to extract the growth 
    # medium label as given in the metadata_table file, again using the strain
    # label which should match between the two files. 

    # assume everything up to '.mz' is the identifier/label of this strain
    strain_name = hdr[:hdr.index('.mz')]
    growth_medium = None 

    # check if the strain label exists in the set of strain mappings the user
    # has provided
    if strain_name not in strains:
        # if this check fails, it could mean a missing strain ID mapping. this isn't
        # fatal but will produce a warning and an entry added to the file 
        # unknown_strains_met.csv in the dataset folder. 
        # 
        # however there could be some identifiers which are not strains and 
        # these should just be ignored. 
        # 
        # if we have a metadata table file, and the parsed strain name does NOT
        # appear in it, this indicates it's not a strain
        if md_table is not None and strain_name not in md_table:
            # we can ignore this completely
            return (None, None, False)

        # if the strain is in the table, then it means the user probably hasn't given us
        # a valid mapping for it, so a warning will be displayed and the name added
        # to the unknown_strains.csv file later on. however if the config file option
        # extended_metadata_table_parsing has been enabled (the ext_metadata_parsing) 
        # parameter in this method), it means we should take the content of the 
        # ATTRIBUTE_Strain column from the table and use that as the strain, adding
        # an extra alias for it too. Additionally the ATTRIBUTE_Medium column content
        # should be recorded as the growth medium for the strain
        if md_table is not None and strain_name in md_table and ext_metadata_parsing:
            growth_medium = md_table[strain_name]['growthmedium']
            strain_col_content = md_table[strain_name]['strain']

            if strain_col_content in strains:
                strain = strains.lookup(strain_col_content)
                strain.add_alias(strain_name)
                # this merges the set of aliases back into the internal
                # lookup table in the StrainCollection
                strains.add(strain)
            else:
                # if this still fails it's an unknown strain 
                logger.warning('Unknown strain identifier: {} (parsed from {})'.format(strain_name, hdr))
                return (strain_name, growth_medium, True)
        else:
            # emit a warning message to indicate that the user needs to add this to 
            # their strain_mappings file
            logger.warning('Unknown strain identifier: {} (parsed from {})'.format(strain_name, hdr))

    return (strain_name, growth_medium, strain_name not in strains)

def load_clusterinfo_old(gnps_format, strains, filename, spec_dict):
    # each line of this file represents a metabolite. 
    # columns representing strain IDs are *ignored* here in favour of parsing
    # .mz(X)ML filenames from either the AllFiles or UniqueFileSources column.
    # both of these list the .mz(X)ML files the molecule was found in (plus the scan 
    # number in the file in the AllFiles case)
    unknown_strains = {}
    with open(filename, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        clu_index_index = headers.index('cluster index')
        if gnps_format == GNPS_FORMAT_OLD_ALLFILES:
            mzxml_index = headers.index('AllFiles')
        elif gnps_format == GNPS_FORMAT_OLD_UNIQUEFILES:
            mzxml_index = headers.index('UniqueFileSources')
        else:
            raise Exception('Unexpected GNPS format {}'.format(gnps_format))

        for line in reader:
            # get the values of the important columns
            clu_index = int(line[clu_index_index]) 
            if gnps_format == GNPS_FORMAT_OLD_UNIQUEFILES:
                mzxmls = line[mzxml_index].split('|')
            else:
                mzxmls = line[mzxml_index].split('###')

            metadata = {'cluster_index': clu_index, 'files': {}}
            seen_files = set()

            for data in mzxmls:
                # TODO ok to ignore scan number if available?
                mzxml = data if gnps_format == GNPS_FORMAT_OLD_UNIQUEFILES else data.split(':')[0]

                # TODO is this OK/sensible?
                if mzxml in seen_files:
                    continue
                seen_files.add(mzxml)

                # add to the list of files for this molecule
                metadata['files'][mzxml] = mzxml

                # should have a strain alias for the mxXML file to lookup here (in theory at least)
                strain = messy_strain_naming_lookup(mzxml, strains)
                if strain is None:
                    # TODO: how to handle this? can't just give up, simply warn?
                    if mzxml not in unknown_strains:
                        logger.warning('Unknown strain: {} for cluster index {}'.format(mzxml, clu_index))
                        unknown_strains[mzxml] = 1
                    else:
                        unknown_strains[mzxml] += 1
                # else:
                #     print('{} ===> {}'.format(mzxml, strain))

                if strain is not None:
                    # TODO this need fixed somehow (missing growth medium info)
                    spec_dict[clu_index].add_strain(strain, None, None, 1)
                
                # update metadata on Spectrum object
                spec_dict[clu_index].metadata.update(metadata)

    if len(unknown_strains) > 0:
        logger.warning('{} unknown strains were detected a total of {} times'.format(len(unknown_strains), sum(unknown_strains.values())))

    return unknown_strains

def parse_metadata_table(filename):
    """Parses the metadata table file from GNPS"""
    if filename is None: 
        return None

    table = {}
    with open(filename, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)

        # check the format is as expected
        if len(headers) < 3 or headers[0] != 'filename':
            # expecting at least 3 columns with the first being 'filename'
            logger.warning('Unrecognised metadata table format in file "{}"'.format(filename))
            return None

        # find the column numbers we're interested in (can't rely on both of these being present)
        # ATTRIBUTE_SampleType, _Strain, _Medium, _Organism
        sampletype_col, medium_col, strain_col = -1, -1, -1
        col_names = {'sampletype': 'ATTRIBUTE_SampleType', 
                    'growthmedium': 'ATTRIBUTE_Medium', 
                    'strain': 'ATTRIBUTE_Strain'}

        # also: ATTRIBUTE_Strain might be useful, there's also ATTRIBUTE_Organism
        for i, hdr in enumerate(headers):
            if hdr == col_names['sampletype']:
                sampletype_col = i
            if hdr == col_names['growthmedium']:
                medium_col = i
            if hdr == col_names['strain']:
                strain_col = i

        for col, name in [(sampletype_col, col_names['sampletype']), 
                (medium_col, col_names['growthmedium']),
                (strain_col, col_names['strain'])]:
            if col == -1:
                logger.warning('No {} column in file "{}"'.format(name, filename))

        for line in reader:
            # only non-BLANK entries
            if sampletype_col != -1 and line[sampletype_col] != 'BLANK':
                # want to strip out the .mz(X)ML extension here to use as the key,
                # then record the available columns (assume SampleType is always
                # available)
                data = {'sampletype': line[sampletype_col], 'growthmedium': None, 'strain': None}
                if medium_col != -1:
                    data['growthmedium'] = line[medium_col]
                if strain_col != -1:
                    data['strain'] = line[strain_col]
                table[RE_MZML_MZXML.sub('', line[0])] = data

    return table

def load_clusterinfo_fbmn(strains, nodes_file, extra_nodes_file, md_table_file, spec_dict, ext_metadata_parsing):
    spec_info = {}

    # parse metadata table if available
    md_table = parse_metadata_table(md_table_file)


    # combine the information in the nodes_file (clusterinfo_summary folder) and
    # the extra_nodes_file (quantification_table folder), indexed by the "cluster index" 
    # and "row ID" fields respectively to link the rows
    with open(nodes_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        ci_index = headers.index('cluster index')

        # take each line and generate a dict storing each header:column value 
        # and insert that dict into the spec_info dict
        for line in reader:
            tmp = {}
            for i, v in enumerate(line):
                tmp[headers[i]] = v
            spec_info[int(line[ci_index])] = tmp

    with open(extra_nodes_file, 'r') as f:
        reader = csv.reader(f, delimiter=',')
        headers = next(reader)
        ri_index = headers.index('row ID')

        # do almost the same thing again but a) match the "cluster index" from the
        # nodes_file to the "row ID" from this file, and update the per-row dict
        # with the extra columns from this file
        for line in reader:
            ri = int(line[ri_index])
            tmp = {}
            for i, v in enumerate(line):
                tmp[headers[i]] = v
            spec_info[ri].update(tmp)

    logger.info('Merged nodes data (new-style), total lines = {}'.format(len(spec_info)))

    unknown_strains = {}

    # TODO: at the moment this iterates over each spectrum in the outer loop, and
    # over the full set of columns in the inner loop. probably makes more sense to 
    # swap that around, or at least to avoid calling parse_mzxml_header repeatedly
    # for the same column headers over and over!

    # for each spectrum 
    for spec_id, spec_data in spec_info.items():
        # look up the Spectrum object with this ID (spec_dict is sourced from the MGF file)
        spectrum = spec_dict[spec_id]
    
        # TODO this should probably be handled by parse_metadata_table 
        if 'ATTRIBUTE_SampleType' in spec_data:
            st = spec_data['ATTRIBUTE_SampleType'].split(',')
            spectrum.metadata['ATTRIBUTE_SampleType'] = st

        # TODO better way of filtering/converting all this stuff down to what's relevant?
        # could search for each strain ID in column title but would be slower?
        # 
        # iterate over the column values for the spectrum
        for k, v in spec_data.items():
            # if the value is a "0", ignore immediately as we only care about nonzero peak intensities
            if v == "0":
                continue

            # ignore any non-"Peak area" columns, want "<strain ID>.mz(X)ML Peak area"
            if k.find(' Peak area') == -1:
                # record the value as an entry in the metadata dict in case it's useful
                # for some other purpose
                spectrum.metadata[k] = v
                continue

            # TODO this will probably need updating for platform data (should already
            # have a set of strain names...)
            #
            # given a "<strain ID>.mz(X)ML Peak area" column value, attempt to lookup the
            # matching Strain object from the set of mappings we have (optionally doing
            # some further stuff involving the metadata table)
            (strain_name, growth_medium, is_unknown) = parse_mzxml_header(k, strains, md_table, ext_metadata_parsing)

            # if the strain name couldn't be recognised at all, ignore it
            if strain_name is None:
                continue 

            # if the name seems valid (appears in the metadata table), record it as unknown
            if is_unknown:
                unknown_strains[k[:k.index('.mz')]] = spec_id
                continue

            # create a new strain object if the intensity value is a float > 0
            v = md_convert(v)
            if strain_name in strains and isinstance(v, float) and v > 0:
                # find the strain object, and add the growth medium + intensity to it. the 
                # growth medium will only be set if extended_metadata_table_parsing is
                # enabled in the config file and the metadata table file contains that info
                strain = strains.lookup(strain_name)
                spectrum.add_strain(strain, strain_name, growth_medium, v)

            # record this as an entry in the metadata dict as well
            spectrum.metadata[k] = v

    return spec_info, unknown_strains

def load_dataset(strains, mgf_file, edges_file, nodes_file, quant_table_file=None, metadata_table_file=None, ext_metadata_parsing=False):
    # common steps to all formats of GNPS data:
    #   - parse the MGF file to create a set of Spectrum objects
    #   - parse the edges file and update the spectra with that data

    # build a set of Spectrum objects by parsing the MGF file
    ms1, ms2, metadata = LoadMGF(name_field='scans').load_spectra([mgf_file])
    logger.info('{} molecules parsed from MGF file'.format(len(ms1)))
    spectra = mols_to_spectra(ms2, metadata)
    # above returns a list, create a dict indexed by spectrum_id to make
    # the rest of the parsing a bit simpler
    spec_dict = {spec.spectrum_id : spec for spec in spectra}

    # add edges info to the spectra
    load_edges(edges_file, spec_dict)

    gnps_format = identify_gnps_format(nodes_file, quant_table_file is not None)
    logger.debug('Nodes_file: {}, quant_table_exists?: {}'.format(nodes_file, quant_table_file is None))
    if gnps_format == GNPS_FORMAT_UNKNOWN:
        raise Exception('Unknown/unsupported GNPS data format')

    # now things depend on the dataset format
    # if we don't have a quantification table, must be older-style dataset (or unknown format)
    if gnps_format != GNPS_FORMAT_NEW_FBMN and quant_table_file is None:
        logger.info('Found older-style GNPS dataset, no quantification table')
        unknown_strains = load_clusterinfo_old(gnps_format, strains, nodes_file, spec_dict)
    else:
        logger.info('quantification table exists, new-style GNPS dataset')
        _, unknown_strains = load_clusterinfo_fbmn(strains, nodes_file, quant_table_file, metadata_table_file, spec_dict, ext_metadata_parsing)

    molfams = make_families(spectra)
    return spec_dict, spectra, molfams, unknown_strains

def make_families(spectra):
    families = []
    family_dict = {}
    family_index = 0
    fams, singles = 0, 0
    for spectrum in spectra:
        family_id = spectrum.family_id
        if family_id == -1: # singleton
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
                new_family.family_id = family_id # preserve original ID here
                family_index += 1
                
                new_family.add_spectrum(spectrum)
                spectrum.family = new_family
                families.append(new_family)
                family_dict[family_id] = new_family
                fams += 1
            else:
                family_dict[family_id].add_spectrum(spectrum)
                spectrum.family = family_dict[family_id]

    logger.debug('make_families: {} molams + {} singletons'.format(fams, singles))
    return families
