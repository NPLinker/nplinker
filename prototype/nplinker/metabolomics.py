import csv

import numpy as np

from .parsers.mgf import LoadMGF
from .strains import StrainCollection
from .annotations import GNPS_KEY
from .utils import sqrt_normalise

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
        # this is a dict indexed by Strain objects (the strains found in this Spectrum), with 
        # the values being dicts of the form {growth_medium: peak intensity} for the parent strain
        self.strains = {}
        self.dataset_strains = None
        self.family = None
        self.random_spectrum = None
        # a dict indexed by filename, or "gnps" 
        self.annotations = {}
        self._losses = None
        self._jcamp = None


    def add_strain(self, strain, growth_medium, peak_intensity):
        if strain in self.strains and growth_medium in self.strains[strain]:
            raise Exception('Clash: {} / {} {}'.format(self, strain, growth_medium))
        
        if strain not in self.strains:
            self.strains[strain] = {}

        self.strains[strain].update({growth_medium: peak_intensity})
        
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

    @property
    def strain_list(self):
        return self.strains.keys()

    @property
    def strain_set(self):
        return set(self.strains.keys())

    def has_strain(self, strain):
        return strain in self.strains

    def get_growth_medium(self, strain):
        if strain not in self.strains:
            return None

        gms = self.strains[strain]
        return list(gms.keys())[0]

    def to_jcamp_str(self, force_refresh=False):
        if self._jcamp is not None and not force_refresh:
            return self._jcamp

        peakdata = '\\n'.join('{}, {}'.format(*p) for p in self.peaks)
        self._jcamp = JCAMP.format(str(self), self.n_peaks, peakdata)
        return self._jcamp

    def add_random(self, strain_prob_dict):
        self.random_spectrum = RandomSpectrum(self, strain_prob_dict)

    def __str__(self):
        return "Spectrum(id={}, spectrum_id={}, strains={})".format(self.id, self.spectrum_id, len(self.strains))

    def __repr__(self):
        return str(self)

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


class RandomSpectrum(object):

    def __init__(self, real_spectrum, strain_prob_dict):
        self.real_spectrum = real_spectrum
        n_strains = 0
        for strain in strain_prob_dict:
            if self.real_spectrum.has_strain(strain):
                n_strains += 1

        self.metadata = set(np.random.choice(list(strain_prob_dict.keys()), n_strains, replace=True, p=list(strain_prob_dict.values())))

    @property
    def strain_list(self):
        return list(self.metadata)

    @property
    def strain_set(self):
        return self.metadata

    def has_strain(self, strain):
        return strain in self.strain_set

def load_spectra(mgf_file):
    ms1, ms2, metadata = LoadMGF(name_field='scans').load_spectra([mgf_file])
    logger.info("load_spectra loaded {} molecules".format(len(ms1)))
    return mols_to_spectra(ms2, metadata)

from .annotations import GNPS_KEY, GNPS_DATA_COLUMNS, create_gnps_annotation

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
            if not GNPS_KEY in new_spectrum.annotations:
                gnps_anno = {k: None for k in GNPS_DATA_COLUMNS}
                gnps_anno['SpectrumID'] = new_spectrum.metadata['spectrumid']
                create_gnps_annotation(new_spectrum, gnps_anno)
        spectra.append(new_spectrum)

    return spectra

def md_convert(val):
    """Try to convert raw metadata values from text to integer, then float"""
    try:
        return int(val)
    except (ValueError, TypeError):
        try:
            return float(val)
        except (ValueError, TypeError):
            if val.lower() == 'n/a':
                return None
    return val

def parse_mzxml_header(hdr, strains, md_table):
    # ignore any non-"Peak area" columns
    # TODO can this be ".mzML" too??
    if hdr.find('.mzXML Peak area') == -1:
        return (None, None)

    # these columns should be named "<strain>_<growth media>.mzXML Peak area"
    #   e.g. KRD178_ISP3.mzXML Peak area
    #   (OR  KRD_178_ISP3.mzXML Peak area)
    # want to extract both strain ID (assumed to be everything up to final underscore)
    # and the growth media
    identifier = hdr[:hdr.index('.')]
    strain_name = None
    growth_medium = '<default>' # TODO what happens if not parsed?
    tokens = identifier.split('_')
    
    # treat the text as a strain ID and carry on)
    if len(tokens) == 1:
        strain_name = identifier
    elif len(tokens) == 2:
        # this might be <strain>_<growthmedia> OR <strain>_<strain_pt2> 
        # (e.g. KR7_178 instead of KRD178)
        if tokens[0] in strains:
            strain_name = tokens[0]
            growth_medium = tokens[1]
        else:
            # check if the original strain is known
            if strain_name not in strains:
                # try removing the underscore
                if ''.join(tokens) in strains:
                    strain_name = ''.join(tokens)
                else:
                    raise Exception('Failed to parse header: {}'.format(hdr))
    elif len(tokens) == 3:
        # this should always(?) be <strain>_<strainpt2>_<growthmedia>
        strain_name = ''.join(tokens[:2])
        growth_medium = tokens[2]
    else:
        raise Exception('Unknown mzXML header format: {}'.format(hdr))

    # found at least one case where these were different cases (m1 vs M1)
    if growth_medium is not None:
        growth_medium = growth_medium.upper()

    # check the final strain_name is valid
    if strain_name not in strains:
        # if this check fails, it could mean a missing strain ID mapping, which should
        # throw an immediate exception. however there could be some identifiers which
        # are not strains and these should just be ignored. 

        # if we have a metadata table file, and the parsed strain name does NOT
        # appear in it, this indicates it's not a strain, so we should ignore it
        if md_table is not None and identifier not in md_table:
            # ignore this completely
            return (None, None)

        # throw an exception in case of unknown strain
        raise Exception('Unknown strain identifier: {} (parsed from {})'.format(strain_name, hdr))

    return (strain_name, growth_medium)

def parse_metadata_table(filename):
    """Parses the metadata table file from GNPS"""
    if filename is None: 
        return None

    table = {}
    with open(filename, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)

        for line in reader:
            # only non-BLANK entries
            if line[1] != 'BLANK':
                table[line[0].replace('.mzXML', '')] = line[1]

    return table

# load metadata from the old-style paired platform TSV format
def _load_metadata_old_platform(strains, nodes_file, spec_dict, growth_media):
    spec_info = {}

    # each line of this file represents a metabolite. the columns representing
    # strain IDs are *ignored* here in favour of parsing the AllFiles column, 
    # which lists the .mzXML files the molecule was found in, plus the scan 
    # number in each file. Also want to get the componentindex column, this is
    # the molecular family (if any, -1 = None?)
    with open(nodes_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        clu_index_index = headers.index('cluster index')
        comp_index_index = headers.index('componentindex')
        allfiles_index = headers.index('AllFiles')

        li = 1
        for line in reader:
            # get the values of the important columns
            clu_index = int(line[clu_index_index]) 
            comp_index = int(line[comp_index_index])
            allfiles = line[allfiles_index].split('###')
            metadata = {'cluster_index': clu_index, 'componentindex': comp_index, 'files': {}}
            seen_files = set()
            for data in allfiles:
                if len(data) < 3:
                    continue

                # TODO: is scan number important for anything?
                # TODO: does it only matter if a molecule appears in a file, not how
                # many times it appears???

                # each entry in the allfiles list will be <mzXML filename>:<scan number>
                mzxml, scan_num = data.split(':')
                scan_num = int(scan_num)
                # TODO re above
                if mzxml in seen_files:
                    continue
                seen_files.add(mzxml)
                metadata['files'][mzxml] = scan_num

                # should have a strain alias for the mxXML file to lookup here
                strain = strains.lookup(mzxml)
                if strain is not None:
                    # TODO growth_medium stuff a bit messy
                    spec_dict[clu_index].add_strain(strain, list(growth_media[strain.id])[0], 1)
            spec_info[clu_index] = metadata
            li += 1

    return spec_info

# load metadata from the old Crusemann or Crusemann-like dataset
def _load_metadata_old(strains, nodes_file, spec_dict):
    spec_info = {}

    with open(nodes_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        ci_index = headers.index('cluster index')

        for line in reader:
            ci = int(line[ci_index])
            metadata = {headers[i]: line[i] for i in range(len(line)) if i != ci_index}
            spec_info[ci] = metadata

    logger.info('Merged nodes data (old-style), total lines = {}'.format(len(spec_info)))

    for spec_id, spec_data in spec_info.items():
        spectrum = spec_dict[spec_id]

        for k, v in spec_data.items():
            # if value is a "0", ignore immediately
            if v == "0":
                continue
            v = md_convert(v)
            # if header indicates this is a strain and it is present in the current
            # spectra, add that info 
            if k in strains:
                strain = strains.lookup(k)
                spectrum.add_strain(strain, None, v)

            spectrum.metadata[k] = v
    return spec_info

def _load_metadata_new(strains, nodes_file, extra_nodes_file, spec_dict, md_table):
    spec_info = {}

    # get a list of the lines in each file, indexed by the "cluster index" and "row ID" fields 
    with open(nodes_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        ci_index = headers.index('cluster index')

        for line in reader:
            tmp = {}
            for i, v in enumerate(line):
                tmp[headers[i]] = v
            spec_info[int(line[ci_index])] = tmp

    with open(extra_nodes_file, 'r') as f:
        reader = csv.reader(f, delimiter=',')
        headers = next(reader)

        ri_index = headers.index('row ID')

        for line in reader:
            ri = int(line[ri_index])
            assert(ri in spec_info)
            tmp = {}
            for i, v in enumerate(line):
                tmp[headers[i]] = v
            spec_info[ri].update(tmp)

    logger.info('Merged nodes data (new-style), total lines = {}'.format(len(spec_info)))

    # for each spectrum 
    for spec_id, spec_data in spec_info.items():
        spectrum = spec_dict[spec_id]
        # TODO better way of filtering/converting all this stuff down to what's relevant?
        # could search for each strain ID in column title but would be slower?
        for k, v in spec_data.items():
            # if the value is a "0", ignore immediately
            if v == "0":
                continue
            (strain_name, growth_medium) = parse_mzxml_header(k, strains, md_table)
            if strain_name is None:
                continue 
            # add a strain object if the value is a float > 0
            v = md_convert(v)
            if strain_name in strains and isinstance(v, float) and v > 0:
                # find the strain object, and add the growth medium + intensity to it
                strain = strains.lookup(strain_name)
                spectrum.add_strain(strain, growth_medium, v)

            spectrum.metadata[k] = v
    return spec_info

def load_metadata(growth_media, strains, nodes_file, extra_nodes_file, metadata_table_file, spectra, spec_dict):
    # this code is complicated by the need to support the Crusemann dataset, which is in a
    # different format to more recent datasets (and any expected future datasets).
    # with the new style datasets, ALL of the following should typically be present:
    #  - nodes_file
    #  - extra_nodes_file
    #  - metadata_table_file 
    # but for crusemann we only have the nodes_file, which has a different set of columns
    # than in the standard case, so it has to be parsed differently 

    if not growth_media:
        # load the rest of the data depending on presence/absence of files
        if extra_nodes_file is None:
            logger.warning('Assuming Crusemann-style dataset!')
            spec_info = _load_metadata_old(strains, nodes_file, spec_dict)
        else:
            # parse metadata table file, if provided
            md_table = parse_metadata_table(metadata_table_file)
            spec_info = _load_metadata_new(strains, nodes_file, extra_nodes_file, spec_dict, md_table)
    else:
        # TODO "new-style" data
        spec_info = _load_metadata_old_platform(strains, nodes_file, spec_dict, growth_media)

    return spectra

def load_edges(edges_file, spectra, spec_dict):

    logger.debug('loading edges file: {} [{} spectra from MGF]'.format(edges_file, len(spectra)))

    with open(edges_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader) # skip headers
        for line in reader:
            spec1_id = int(line[0])
            spec2_id = int(line[1])
            cosine = float(line[4])
            family = int(line[6])

            spec1 = spec_dict[spec1_id]
            spec2 = spec_dict[spec2_id]

            if family != -1: # singletons
                spec1.family = family
                spec2.family = family

                spec1.edges.append((spec2.id, spec2.spectrum_id, cosine))
                spec2.edges.append((spec1.id, spec1.spectrum_id, cosine))
            else:
                spec1.family = family

    return spectra

class MolecularFamily(object):
    def __init__(self, family_id):
        self.id = -1 
        self.family_id = family_id
        self.spectra = []
        self.family = None
        self.random_molecular_family = RandomMolecularFamily(self)

    def has_strain(self, strain):
        for spectrum in self.spectra:
            if spectrum.has_strain(strain):
                return True

        return False

    def add_spectrum(self, spectrum):
        self.spectra.append(spectrum)

    def __str__(self):
        return 'MolFam(family_id={}, spectra={})'.format(self.family_id, len(self.spectra))

class RandomMolecularFamily(object):
    def __init__(self, molecular_family):
        self.molecular_family = molecular_family

    def has_strain(self, strain):
        for spectrum in self.molecular_family.spectra:
            if hasattr(spectrum, 'random_spectrum'):
                if spectrum.random_spectrum.has_strain(strain):
                    return True
            else:
                print("Spectrum objects need random spectra attached for this functionality")

        return False

class SingletonFamily(MolecularFamily):
    def __init__(self):
        super(SingletonFamily, self).__init__(-1)

    def __str__(self):
        return "Singleton molecular family"

# TODO this should be better integrated with rest of the loading process
# so that spec.family is never set to a primitive then replaced by an
# object only if this method is called
# (this is now called every time at least)
def make_families(spectra):
    families = []
    family_dict = {}
    family_index = 0
    fams, singles = 0, 0
    for spectrum in spectra:
        family_id = spectrum.family
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

def read_aa_losses(filename):
    """
    Read AA losses from data file. (assume fixed structure...)
    """
    aa_losses = {}
    with open(filename, 'r') as f:
        reader = csv.reader(f, delimiter=',')
        next(reader) # skip headers
        for line in reader:
            if len(line) == 0:
                continue
            aa_id = line[1]
            aa_mono = float(line[4])
            aa_avg = float(line[5])
            aa_losses[aa_id.lower()] = (aa_mono, aa_avg)

    return aa_losses
