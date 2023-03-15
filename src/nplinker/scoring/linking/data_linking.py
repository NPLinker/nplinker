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

# Methods to find correlations between spectra/molecular families and
# gene clusters/families (BGCs/GCFs)
#
# (still at very much protoype/exploration stage)
#
# Naming:
# M_*   stands for a matrix format
# map_* stands for a simple mapping lookup table
# spec  stands for spectrum
# fam   stands for molecular family

from collections import Counter
from typing import Sequence
# import packages
import numpy as np
import pandas as pd

from nplinker.metabolomics.molecular_family import MolecularFamily
from nplinker.genomics.gcf import GCF
from nplinker.metabolomics.spectrum import Spectrum
from .data_linking_functions import calc_correlation_matrix


SCORING_METHODS = ['metcalf', 'likescore', 'hg']

from nplinker.logconfig import LogConfig


logger = LogConfig.getLogger(__name__)


class DataLinks():
    """
    DataLinks collects and structures co-occurence data
    1) Co-occurences of spectra, families, and GCFs with respect to strains
    2) Mappings: Lookup-tables that link different ids and categories
    3) Correlation matrices that show how often spectra/families and GCFs co-occur
    """

    def __init__(self):
        # matrices that store co-occurences with respect to strains
        # values = 1 where gcf/spec/fam occur in strain
        # values = 0 where gcf/spec/fam do not occur in strain
        self.M_gcf_strain = []
        self.M_spec_strain = []
        self.M_fam_strain = []

        # mappings (lookup lists to map between different ids and categories
        self.mapping_spec = pd.DataFrame()
        self.mapping_gcf = pd.DataFrame()
        self.mapping_fam = pd.DataFrame()  # labels for strain-family matrix
        self.mapping_strain = pd.DataFrame()
        self.family_members = []

        # correlation matrices for spectra <-> GCFs
        self.M_spec_gcf = [
        ]  # = int: Number of strains where spec_x and gcf_y co_occure
        self.M_spec_notgcf = [
        ]  # = int: Number of strains where spec_x and NOT-gcf_y co_occure
        self.M_notspec_gcf = [
        ]  # = int: Number of strains where NOT-spec_x and gcf_y co_occure
        # and the same for mol.families <-> GCFs
        self.M_fam_gcf = []
        self.M_fam_notgcf = []
        self.M_notfam_gcf = []

    def get_spec_pos(self, spec_id):
        # get the position in the arrays of a spectrum
        row = self.mapping_spec.loc[self.mapping_spec['original spec-id'] ==
                                    float(spec_id)]
        return int(row.iloc[0]['spec-id'])

    def get_gcf_pos(self, gcf_id):
        # TODO: fix this so the original ID is present in case of re-ordering
        pass

    def load_data(self, spectra, gcf_list, strain_list, molfams):
        # load data from spectra, GCFs, and strains
        logger.debug("Create mappings between spectra, gcfs, and strains.")
        self.collect_mappings_spec(spectra)
        # self.collect_mappings_spec_v2(molfams)
        self.collect_mappings_gcf(gcf_list)
        logger.debug(
            "Create co-occurence matrices: spectra<->strains + and gcfs<->strains."
        )
        self.matrix_strain_gcf(gcf_list, strain_list)
        self.matrix_strain_spec(spectra, strain_list)

    def find_correlations(self, include_singletons=False):
        # collect correlations/ co-occurences
        logger.debug("Create correlation matrices: spectra<->gcfs.")
        self.correlation_matrices(type='spec-gcf')
        logger.debug("Create correlation matrices: mol-families<->gcfs.")
        self.data_family_mapping(include_singletons=include_singletons)
        self.correlation_matrices(type='fam-gcf')

    def collect_mappings_spec(self, obj: Sequence[Spectrum]|Sequence[MolecularFamily]):
        if isinstance(obj[0], Spectrum):
            mapping_spec = self._collect_mappings_from_spectra(obj)
        elif isinstance(obj[0], MolecularFamily):
            mapping_spec = self._collect_mappings_from_molecular_families(obj)

        # extend mapping tables:
        self.mapping_spec["spec-id"] = mapping_spec[:, 0]
        self.mapping_spec["original spec-id"] = mapping_spec[:, 1]
        self.mapping_spec["fam-id"] = mapping_spec[:, 2]

    def _collect_mappings_from_spectra(self, spectra) -> np.ndarray[np.float64]:
        # Collect most import mapping tables from input data
        mapping_spec = np.zeros((len(spectra), 3))
        mapping_spec[:, 0] = np.arange(0, len(spectra))

        for i, spectrum in enumerate(spectra):
            mapping_spec[i, 1] = spectrum.id
            mapping_spec[i, 2] = spectrum.family.family_id

        return mapping_spec

    def _collect_mappings_from_molecular_families(self, molfams: Sequence[MolecularFamily]) -> np.ndarray[np.float64]:
        num_spectra = sum(len(x.spectra_ids) for x in molfams)
        mapping_spec = np.zeros((num_spectra, 3))
        mapping_spec[:, 0] = np.arange(0, num_spectra)

        inverted_mappings = {}
        for item in molfams:
            for spectrum_id in item.spectra_ids:
                inverted_mappings[spectrum_id] = item.family_id

        for i, key in enumerate(inverted_mappings):
            mapping_spec[i, 1] = key
            mapping_spec[i, 2] = inverted_mappings[key]

        return mapping_spec

    def collect_mappings_gcf(self, gcf_list):
        """
        Find classes of gene cluster (nrps, pksi etc.)
        collect most likely class (most occuring name, preferentially not "Others")
        additional score shows fraction of chosen class among all given ones
        """

        # TODO: not only collect bigclass types but also product predictions
        bigscape_bestguess = []
        for i, gcf in enumerate(gcf_list):
            #bigscape_class = []
            #for i, bgc in enumerate(gcf_list[i].bgcs):
            #    # bigscape_class.append(gcf_list[i].bgc_list[m].bigscape_class)
            #    bigscape_class.append(bgc.bigscape_class)
            #    class_counter = Counter(bigscape_class)

            # TODO: this might need properly rewritten due to changes in the way GCF/BGC
            # objects store bigscape class information and handle hybrid BGCs (see genomics.py). the
            # original version i've left above iterates over every BGC in the current GCF
            # and extracts its .bigscape_class attribute, but now each BGC can have multiple
            # classes if it happens to be a hybrid and i'm not sure the version below
            # still makes sense.
            #
            # doesn't seem to be very important for calculating the metcalf scores though so not urgent for now.
            bigscape_class = []
            for bgc in gcf.bgcs:
                bigscape_class.extend(bgc.bigscape_classes)
            class_counter = Counter(bigscape_class)

            # try not to select "Others":
            if class_counter.most_common(1)[0][0] is None:
                bigscape_bestguess.append(["Others", 0])
            elif class_counter.most_common(
                    1)[0][0] == "Others" and class_counter.most_common(
                        1)[0][1] < len(bigscape_class):
                if class_counter.most_common(2)[1][0] is None:
                    bigscape_bestguess.append([
                        class_counter.most_common(1)[0][0],
                        class_counter.most_common(1)[0][1] /
                        len(bigscape_class)
                    ])
                else:
                    bigscape_bestguess.append([
                        class_counter.most_common(2)[1][0],
                        class_counter.most_common(2)[1][1] /
                        len(bigscape_class)
                    ])
            else:
                bigscape_bestguess.append([
                    class_counter.most_common(1)[0][0],
                    class_counter.most_common(1)[0][1] / len(bigscape_class)
                ])

        # extend mapping tables:
        self.mapping_gcf["gcf-id"] = np.arange(0, len(bigscape_bestguess))
        bigscape_guess, bigscape_guessscore = zip(*bigscape_bestguess)
        self.mapping_gcf["bgc class"] = bigscape_guess
        self.mapping_gcf["bgc class score"] = bigscape_guessscore

    def matrix_strain_gcf(self, gcf_list, strain_list):
        # Collect co-ocurences in M_spec_strain matrix
        M_gcf_strain = np.zeros((len(gcf_list), len(strain_list)))

        for i, strain in enumerate(strain_list):
            for m, gcf in enumerate(gcf_list):
                if gcf.has_strain(strain):
                    M_gcf_strain[m, i] = 1

        self.M_gcf_strain = M_gcf_strain
        # extend mapping tables:
        self.mapping_gcf["no of strains"] = np.sum(self.M_gcf_strain, axis=1)
        self.mapping_strain["no of gcfs"] = np.sum(self.M_gcf_strain, axis=0)

    def matrix_strain_spec(self, spectra, strain_list):
        # Collect co-ocurences in M_strains_specs matrix

        M_spec_strain = np.zeros((len(spectra), len(strain_list)))
        for i, spectrum in enumerate(spectra):
            for j, s in enumerate(strain_list):
                if spectrum.has_strain(s):
                    M_spec_strain[i, j] = 1
        self.M_spec_strain = M_spec_strain

        # extend mapping tables:
        self.mapping_spec["no of strains"] = np.sum(self.M_spec_strain, axis=1)
        self.mapping_strain["no of spectra"] = np.sum(self.M_spec_strain,
                                                      axis=0)
        self.mapping_strain["strain name"] = [str(s) for s in strain_list]

    def data_family_mapping(self, include_singletons=False):
        # Create M_fam_strain matrix that gives co-occurences between mol. families and strains
        # matrix dimensions are: number of families  x  number of strains

        family_ids = np.unique(
            self.mapping_spec["fam-id"])  # get unique family ids

        # if singletons are included, check if there are a non-zero number of them (singleton
        # families all have -1 as a family ID number)
        if include_singletons and np.where(
                self.mapping_spec["fam-id"] == -1)[0].shape[0] > 0:
            # in this case the number of unique families is the number of singletons
            # plus the number of normal families. the "-1" is (I think) to account for
            # the single "-1" entry that will be present in "family_ids".
            num_of_singletons = np.where(
                self.mapping_spec["fam-id"] == -1)[0].shape[0]
            num_unique_fams = num_of_singletons + len(family_ids) - 1
        else:
            # if no singletons included or present in the dataset, just take the number
            # of regular molfams instead
            num_of_singletons = 0
            num_unique_fams = len(family_ids)

        M_fam_strain = np.zeros((num_unique_fams, self.M_spec_strain.shape[1]))
        strain_fam_labels = []
        strain_fam_index = []

        if num_of_singletons > 0:  # if singletons exist + included
            M_fam_strain[(
                num_unique_fams -
                num_of_singletons):, :] = self.M_spec_strain[np.where(
                    self.mapping_spec["fam-id"][:, 0] == -1)[0], :]

        # go through families (except singletons) and collect member strain occurences
        self.family_members = []
        for i, fam_id in enumerate(
                family_ids[np.where(family_ids != -1)].astype(int)):
            family_members = np.where(
                np.array(self.mapping_spec["fam-id"]) == fam_id)
            self.family_members.append(family_members)
            M_fam_strain[i, :] = np.sum(self.M_spec_strain[family_members, :],
                                        axis=1)
            strain_fam_labels.append(fam_id)
            strain_fam_index.append(i)

        add_singleton_entries = -1 in family_ids
        # TODO: i think this breaks stuff below due to mismatches in the number of rows
        # in the dataframes and matrices if there are no -1 family ids.
        # discovered when trying to write some code to test scoring. is this ever
        # likely to happen with a real dataset??
        if add_singleton_entries:
            strain_fam_labels.append([-1] * num_of_singletons)
            strain_fam_index.append(i + 1)

        # only looking for co-occurence, hence only 1 or 0
        M_fam_strain[M_fam_strain > 1] = 1

        self.M_fam_strain = M_fam_strain
        # extend mapping table:
        self.mapping_fam["family id"] = strain_fam_index
        self.mapping_fam["original family id"] = strain_fam_labels
        self.mapping_fam["no of strains"] = np.sum(self.M_fam_strain, axis=1)
        num_members = [x[0].shape[0] for x in self.family_members]
        # see above
        if add_singleton_entries:
            num_members.append(num_of_singletons)
        self.mapping_fam["no of members"] = num_members
        return self.family_members

    def common_strains(self, objects_a, objects_b, filter_no_shared=False):
        """
        Obtain the set of common strains between all pairs from the lists objects_a
        and objects_b.

        The two parameters can be either lists or single instances of the 3 supported
        object types (GCF, Spectrum, MolecularFamily). It's possible to use a single
        object together with a list as well.

        Returns a dict indexed by tuples of (Spectrum/MolecularFamily, GCF), where
        the values are lists of strain indices which appear in both objects, which
        can then be looked up in NPLinker.strains.
        """

        # TODO make this work for BGCs too?

        is_list_a = isinstance(objects_a, list)
        is_list_b = isinstance(objects_b, list)

        type_a = type(objects_a[0]) if is_list_a else type(objects_a)
        type_b = type(objects_b[0]) if is_list_b else type(objects_b)

        if type_a == type_b:
            raise Exception('Must supply objects with different types!')

        # to keep things slightly simpler, ensure the GCFs are always "b"
        if type_a == GCF:
            type_a, type_b = type_b, type_a
            is_list_a, is_list_b = is_list_b, is_list_a
            objects_a, objects_b = objects_b, objects_a

        if not is_list_a:
            objects_a = [objects_a]
        if not is_list_b:
            objects_b = [objects_b]

        # retrieve object IDs
        # TODO: issue #103 stop using gcf.id, but note that the ids_b should be
        # a list of int
        ids_b = [gcf.id for gcf in objects_b]
        # these might be MolFams or Spectra, either way they'll have a .id attribute
        ids_a = [obj.id for obj in objects_a]

        data_a = self.M_spec_strain if type_a == Spectrum else self.M_fam_strain
        data_b = self.M_gcf_strain

        results = {}
        for a, obj_a in enumerate(objects_a):
            for b, obj_b in enumerate(objects_b):
                # just AND both arrays and extract the indices with positive results
                result = np.where(
                    np.logical_and(data_a[ids_a[a]], data_b[ids_b[b]]))[0]
                # if we want to exclude results with no shared strains
                if (filter_no_shared
                        and len(result) > 0) or not filter_no_shared:
                    results[(obj_a, obj_b)] = result

        return results

    def correlation_matrices(self, type='spec-gcf'):
        """
        Collect co-occurrences accros strains:
        IF type='spec-gcf':
            number of co-occurences of spectra and GCFS
            --> Output: M_spec_gcf matrix
        IF type='fam-gcf':
            number of co-occurences of mol.families and GCFS
            --> Output: M_fam_gcf matrix
        """

        # Make selection for scenario spec<->gcf or fam<->gcf
        if type == 'spec-gcf':
            M_type1_strain = self.M_spec_strain
        elif type == 'fam-gcf':
            M_type1_strain = self.M_fam_strain
        elif type == 'spec-bgc' or type == 'fam-bgc':
            raise Exception("Given types are not yet supported... ")
        else:
            raise Exception(
                "Wrong correlation 'type' given. Must be one of 'spec-gcf', 'fam-gcf', ..."
            )

        logger.debug(
            f"Calculating correlation matrices of type: {type}")

        # Calculate correlation matrix from co-occurence matrices
        M_type1_gcf, M_type1_notgcf, M_nottype1_gcf, M_nottype1_notgcf = calc_correlation_matrix(
            M_type1_strain, self.M_gcf_strain)

        # return results:
        if type == 'spec-gcf':
            self.M_spec_gcf = M_type1_gcf
            self.M_spec_notgcf = M_type1_notgcf
            self.M_notspec_gcf = M_nottype1_gcf
            self.M_notspec_notgcf = M_nottype1_notgcf
        elif type == 'fam-gcf':
            self.M_fam_gcf = M_type1_gcf
            self.M_fam_notgcf = M_type1_notgcf
            self.M_notfam_gcf = M_nottype1_gcf
            self.M_notfam_notgcf = M_nottype1_notgcf
        else:
            raise Exception("No correct correlation matrix was created.")

    # class data_links OUTPUT functions
    # TODO add output functions (e.g. to search for mappings of individual specs, gcfs etc.)
