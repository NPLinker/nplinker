import csv
import numpy as np

from parsers import LoadMGF

from logconfig import LogConfig
logger = LogConfig.getLogger(__file__)

class Spectrum(object):
    
    METADATA_BLACKLIST = set(['AllOrganisms', 'LibraryID', 'RTStdErr', 'RTMean', 'AllGroups', 'DefaultGroups',
                            'precursor mass', 'parent mass', 'ProteoSAFeClusterLink', 'precursor intensity', 'sum(precursor intensity)',
                              'precursormass', 'parentmass', 'singlechargeprecursormass', 'cluster index', 'number of spectra', 
                              'UniqueFileSourcesCount', 'EvenOdd', 'charge', 'precursor charge'])

    def __init__(self, id, peaks, spectrum_id, precursor_mz, parent_mz=None, rt=None):
        self.id = id
        self.peaks = sorted(peaks, key=lambda x: x[0]) # ensure sorted by mz
        self.n_peaks = len(self.peaks)
        self.max_ms2_intensity = max([intensity for mz, intensity in self.peaks])
        self.total_ms2_intensity = sum([intensity for mz, intensity in self.peaks])
        self.spectrum_id = str(spectrum_id)
        self.rt = rt
        self.precursor_mz = precursor_mz
        self.parent_mz = parent_mz
        self.metadata = {}
        self.edges = []
        self.family = None
        self.random_spectrum = None
        self.annotations = []

        self._losses = None

    def annotation_from_metadata(self):
        annotation = self.get_metadata_value('LibraryID')
        if annotation:
            self.annotations.append((annotation, 'gnps'))

    def get_metadata_value(self, key):
        val = self.metadata.get(key, None)
        return val

    @property
    def strain_list(self):
        return [k for k in self.metadata.keys() if self.has_strain(k)]

    @property
    def strain_set(self):
        return set([k for k in self.metadata.keys() if self.has_strain(k)])

    def has_strain(self, strain):
        if strain in Spectrum.METADATA_BLACKLIST:
            return False

        strain_val = self.metadata.get(strain, 0)
        #  TODO this can throw an exception if a value is a non-integer or None...
        try:
            if strain_val > 0:
                return True
        except TypeError:
            # logger.debug('Warning: has_strain({}) failed because metadata.get returned "{}" (type({}))'.format(strain, strain_val, type(strain_val)))
            pass
        return False

    # def print_spectrum(self):
    #     print()
    #     print(self.file_name,self.spectrum_id)
    #     for i,(mz,intensity) in enumerate(self.peaks):
    #         print(i,mz,intensity,self.normalised_peaks[i][1])

    # def plot(self,xlim = None,**kwargs):
    #     plot_spectrum(self.peaks,xlim=xlim,title = "{} {} (m/z= {})".format(self.file_name,self.scan_number,self.parent_mz),**kwargs)


    def add_random(self, strain_prob_dict):
        self.random_spectrum = RandomSpectrum(self, strain_prob_dict)

    def __str__(self):
        # return "Spectrum {} with {} peaks, max_ms2_intensity {}".format(self.spectrum_id,self.n_peaks,self.max_ms2_intensity)
        return "Spectrum(id={}, spectrum_id={})".format(self.id, self.spectrum_id)

    def __repr__(self):
        return str(self)

    def __cmp__(self, other):
        if self.parent_mz >= other.parent_mz:
            return 1
        else:
            return -1

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
    print("Loaded {} molecules".format(len(ms1)))
    return mols_to_spectra(ms2, metadata)

def mols_to_spectra(ms2, metadata):
    ms2_dict = {}
    for m in ms2:
        if not m[3] in ms2_dict:
            ms2_dict[m[3]] = []
        ms2_dict[m[3]].append((m[0], m[2]))
    
    spectra = []
    for i, m in enumerate(ms2_dict.keys()):
        new_spectrum = Spectrum(i, ms2_dict[m], m.name, metadata[m.name]['precursormass'])
        new_spectrum.metadata = metadata[m.name]
        spectra.append(new_spectrum)

    return spectra

def load_additional_annotations(spectra,annotation_file,id_field,annotation_field):
    with open(annotation_file,'r') as f:
        reader = csv.reader(f,delimiter = '\t')
        heads = next(reader)
        source_id_idx = heads.index(id_field)
        annotation_idx = heads.index(annotation_field)
        new_annotations = {}
        for line in reader:
            new_annotations[line[source_id_idx]] = (line[annotation_idx],annotation_file)
    for s in spectra:
        found_comp = set()
        if s.spectrum_id in new_annotations:
            compound = new_annotations[s.spectrum_id][0]
            if not compound in found_comp:
                s.annotations.append(new_annotations[s.spectrum_id])
                found_comp.add(new_annotations[s.spectrum_id][0])
            

def load_metadata(spectra, metadata_file):
    # make a dictionary mapping spectum ids to spectrum
    spec_dict = {}
    for spec in spectra:
        spec_dict[spec.spectrum_id] = spec

    # to collect all strains from this source for later use
    strains = set()

    with open(metadata_file, 'rU') as f:
        reader = csv.reader(f, delimiter='\t')
        heads = next(reader)
        for line in reader:
            if len(line) == 0:
                continue
            spectrum = spec_dict[line[0]]
            for i, value in enumerate(line):
                key = heads[i]
                try:
                    value = int(value)
                except ValueError:
                    if value == 'N/A':
                        value = None
                spectrum.metadata[key] = value
                strains.add(key)
            spectrum.annotation_from_metadata()

    return strains

def load_edges(spectra, edge_file):
    spec_dict = {}
    for spec in spectra:
        spec_dict[spec.spectrum_id] = spec

    with open(edge_file, 'rU') as f:
        reader = csv.reader(f, delimiter='\t')
        heads = next(reader)
        for line in reader:
            spec1_id = line[0]
            spec2_id = line[1]
            cosine = float(line[4])
            family = line[6]

            spec1 = spec_dict[spec1_id]
            spec2 = spec_dict[spec2_id]

            if family != '-1': # singletons
                spec1.family = family
                spec2.family = family

                spec1.edges.append((spec2.spectrum_id, cosine))
                spec2.edges.append((spec1.spectrum_id, cosine))
            else:
                spec1.family = family


class MolecularFamily(object):
    def __init__(self, family_id):
        self.id = -1 
        self.family_id = family_id
        self.spectra = []
        self.random_molecular_family = RandomMolecularFamily(self)

    def has_strain(self, strain):
        for spectrum in self.spectra:
            if spectrum.has_strain(strain):
                return True

        return False

    def add_spectrum(self, spectrum):
        self.spectra.append(spectrum)

    def __str__(self):
        return "Molecular family with {} spectra".format(len(self.spectra))

class RandomMolecularFamily(object):
    def __init__(self, molecular_family):
        self.molecular_family = molecular_family

    def has_strain(self, strain):
        for spectrum in self.molecular_family.spectra:
            try:
                if spectrum.random_spectrum.has_strain(strain):
                    return True
            except:
                print("Spectrum objects need random spectra attached for this functionality")

        return False

class SingletonFamily(MolecularFamily):
    def __init__(self):
        super(SingletonFamily, self).__init__('-1')

    def __str__(self):
        return "Singleton molecular family"

def make_families(spectra):
    families = []
    family_dict = {}
    family_index = 0
    for spectrum in spectra:
        family_id = spectrum.family
        if family_id == '-1': # singleton
            new_family = SingletonFamily()
            new_family.id = family_index
            family_index += 1

            new_family.add_spectrum(spectrum)
            spectrum.family = new_family
            families.append(new_family)
        else:
            if family_id not in family_dict:
                new_family = MolecularFamily(family_id)
                new_family.id = family_index
                family_index += 1
                
                new_family.add_spectrum(spectrum)
                spectrum.family = new_family
                families.append(new_family)
                family_dict[family_id] = new_family
            else:
                family_dict[family_id].add_spectrum(spectrum)
                spectrum.family = family_dict[family_id]


    return families

def read_aa_losses(filename):
    """
    Read AA losses from data file. (assume fixed structure...)
    """
    aa_losses = {}
    with open(filename, 'rU') as f:
        reader = csv.reader(f, delimiter=',')
        header = next(reader)
        for line in reader:
            if len(line) == 0:
                continue
            aa_id = line[1]
            aa_mono = float(line[4])
            aa_avg = float(line[5])
            aa_losses[aa_id.lower()] = (aa_mono, aa_avg)

    return aa_losses
