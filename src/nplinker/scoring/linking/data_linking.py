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

from __future__ import annotations
from collections import Counter
from typing import Sequence, TYPE_CHECKING
import numpy as np
import pandas as pd
from nplinker.genomics.gcf import GCF
from nplinker.logconfig import LogConfig
from nplinker.metabolomics.molecular_family import MolecularFamily
from nplinker.metabolomics.singleton_family import SingletonFamily
from nplinker.metabolomics.spectrum import Spectrum
from .data_linking_functions import calc_correlation_matrix


if TYPE_CHECKING:
    from nplinker.strain_collection import StrainCollection

logger = LogConfig.getLogger(__name__)


class DataLinks():
    """
    DataLinks collects and structures co-occurence data
    1) Co-occurences of spectra, families, and GCFs with respect to strains
    2) Mappings: Lookup-tables that link different ids and categories
    3) Correlation matrices that show how often spectra/families and GCFs co-occur
    """

    def __init__(self):
        # DataFrame to store occurence of objects with respect to strains
        # values = 1 where gcf/spec/fam occur in strain, 0 otherwise
        self.gcf_strain_occurrence = pd.DataFrame()
        self.spec_strain_occurrence = pd.DataFrame()
        self.mf_strain_occurrence = pd.DataFrame()

        # mappings (lookup lists to map between different ids and categories
        self.mapping_spec = pd.DataFrame()
        self.mapping_gcf = pd.DataFrame()
        self.mapping_fam = pd.DataFrame()  # labels for strain-family matrix
        self.mapping_strain = pd.DataFrame()

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

    def load_data(self, spectra: Sequence[Spectrum], gcfs: Sequence[GCF],
                  strains: StrainCollection,
                  molfams: Sequence[MolecularFamily]):
        # load data from spectra, GCFs, and strains
        logger.debug("Create mappings between spectra, gcfs, and strains.")
        self.collect_mappings_spec(spectra)
        self.collect_mappings_gcf(gcfs)
        logger.debug(
            "Create co-occurence matrices: spectra<->strains, gcfs<->strains and mfs<->strains."
        )
        self._get_gcf_strain_occurrence(gcfs, strains)
        self._get_spec_strain_occurrence(spectra, strains)
        self._get_mf_strain_occurrence(molfams, strains)

        self._get_mappings_from_occurrence()

    def _get_mappings_from_occurrence(self):
        self.mapping_gcf["no of strains"] = np.sum(self.gcf_strain_occurrence,
                                                   axis=1)
        self.mapping_spec["no of strains"] = np.sum(
            self.spec_strain_occurrence, axis=1)
        self.mapping_strain["no of spectra"] = np.sum(
            self.spec_strain_occurrence, axis=0)
        self.mapping_fam["no of strains"] = np.sum(self.mf_strain_occurrence,
                                                   axis=1)

    def find_correlations(self):
        # collect correlations/ co-occurences
        logger.debug("Create correlation matrices: spectra<->gcfs.")
        self.correlation_matrices(type='spec-gcf')
        logger.debug("Create correlation matrices: mol-families<->gcfs.")
        self.correlation_matrices(type='fam-gcf')

    def collect_mappings_spec(self, obj: Sequence[Spectrum]
                              | Sequence[MolecularFamily]):
        if isinstance(obj[0], Spectrum):
            mapping_spec = self._collect_mappings_from_spectra(obj)
        elif isinstance(obj[0], MolecularFamily):
            mapping_spec = self._collect_mappings_from_molecular_families(obj)

        # extend mapping tables:
        # TODO: why do we need the mappings???
        # "spec-id" is defined as the index of the spectrum in the input data
        self.mapping_spec["fam-id"] = mapping_spec[:, 2]

    def _collect_mappings_from_spectra(self,
                                       spectra) -> np.ndarray[np.float64]:
        # Collect most import mapping tables from input data
        mapping_spec = np.zeros((len(spectra), 3))
        mapping_spec[:, 0] = np.arange(0, len(spectra))

        for i, spectrum in enumerate(spectra):
            mapping_spec[i, 1] = spectrum.id
            mapping_spec[i, 2] = spectrum.family.family_id

        return mapping_spec

    def _collect_mappings_from_molecular_families(
            self,
            molfams: Sequence[MolecularFamily]) -> np.ndarray[np.float64]:
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
        bigscape_guess, bigscape_guessscore = zip(*bigscape_bestguess)
        self.mapping_gcf["bgc class"] = bigscape_guess

    def _get_gcf_strain_occurrence(self, gcfs: Sequence[GCF],
                                   strains: StrainCollection) -> None:
        """Get the occurence of strains in gcfs.

        The occurence is a DataFrame with gcfs as rows and strains as columns,
        where index is `gcf.gcf_id` and column name is `strain.id`. The values
        are 1 if the gcf contains the strain and 0 otherwise.
        """
        df_gcf_strain = pd.DataFrame(np.zeros((len(gcfs), len(strains))),
                                     index=[gcf.gcf_id for gcf in gcfs],
                                     columns=[strain.id for strain in strains])
        for gcf in gcfs:
            for strain in strains:
                if gcf.has_strain(strain):
                    df_gcf_strain.loc[gcf.gcf_id, strain.id] = 1
        self.gcf_strain_occurrence = df_gcf_strain

    def _get_spec_strain_occurrence(self, spectra: Sequence[Spectrum],
                                    strains: StrainCollection) -> None:
        """Get the occurence of strains in spectra.

        The occurence is a DataFrame with spectra as rows and strains as columns,
        where index is `spectrum.spectrum_id` and column name is `strain.id`.
        The values are 1 if the spectrum contains the strain and 0 otherwise.
        """
        df_spec_strain = pd.DataFrame(
            np.zeros((len(spectra), len(strains))),
            index=[spectrum.spectrum_id for spectrum in spectra],
            columns=[strain.id for strain in strains])
        for spectrum in spectra:
            for strain in strains:
                if spectrum.has_strain(strain):
                    df_spec_strain.loc[spectrum.spectrum_id, strain.id] = 1
        self.spec_strain_occurrence = df_spec_strain

    def _get_mf_strain_occurrence(self, mfs: Sequence[MolecularFamily],
                                  strains: StrainCollection) -> None:
        """Get the occurence of strains in molecular families.

        The occurence is a DataFrame with molecular families as rows and
        strains as columns, where index is `mf.family_id` and column name is
        `strain.id`. The values are 1 if the molecular family contains the
        strain and 0 otherwise.

        Note that SingletonFamily objects are excluded from given `mfs`.
        """
        # remove SingletonFamily objects
        mfs = [mf for mf in mfs if not isinstance(mf, SingletonFamily)]

        df_mf_strain = pd.DataFrame(np.zeros((len(mfs), len(strains))),
                                    index=[mf.family_id for mf in mfs],
                                    columns=[strain.id for strain in strains])
        for mf in mfs:
            for strain in strains:
                if mf.has_strain(strain):
                    df_mf_strain.loc[mf.family_id, strain.id] = 1
        self.mf_strain_occurrence = df_mf_strain

    def get_common_strains(
        self,
        spectra_or_molfams: Sequence[Spectrum] | Sequence[MolecularFamily],
        gcfs: Sequence[GCF],
        filter_no_shared: bool = False
    ) -> dict[tuple[Spectrum | MolecularFamily, GCF], list[str]]:
        """Get common strains between given spectra/molecular families and GCFs.

        Note that SingletonFamily objects are excluded from given `spectra_or_molfams`.

        Args:
            spectra_or_molfams(Sequence[Spectrum] | Sequence[MolecularFamily]):
                A list of Spectrum or MolecularFamily objects.
            gcfs(Sequence[GCF]): A list of GCF objects.
            filter_no_shared(bool): If True, return only the pair of
                spectrum/molecular family and GCF that have common strains;
                otherwise, return all pairs no matter they have common strains
                or not.

        Returns:
            dict: A dict where the keys are tuples of (Spectrum/MolecularFamily, GCF)
            and values are a list of strain ids that appear in both objects.
        """
        if not isinstance(spectra_or_molfams[0], (Spectrum, MolecularFamily)):
            raise ValueError(
                'Must provide Spectra or MolecularFamilies as the first argument!'
            )
        if not isinstance(gcfs[0], GCF):
            raise ValueError('Must provide GCFs as the second argument!')

        # Assume that 3 occurrence dataframes have same df.columns (strain ids)
        strain_ids = self.gcf_strain_occurrence.columns
        results = {}
        for obj in spectra_or_molfams:
            if isinstance(obj, SingletonFamily):
                continue
            for gcf in gcfs:
                if isinstance(obj, Spectrum):
                    shared_strains = strain_ids[np.logical_and(
                        self.spec_strain_occurrence.loc[obj.spectrum_id],
                        self.gcf_strain_occurrence.loc[gcf.gcf_id])]
                else:
                    shared_strains = strain_ids[np.logical_and(
                        self.mf_strain_occurrence.loc[obj.family_id],
                        self.gcf_strain_occurrence.loc[gcf.gcf_id])]
                if filter_no_shared and len(shared_strains) == 0:
                    continue
                results[(obj, gcf)] = shared_strains.to_list()
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
            M_type1_strain = self.spec_strain_occurrence
        elif type == 'fam-gcf':
            M_type1_strain = self.mf_strain_occurrence
        elif type == 'spec-bgc' or type == 'fam-bgc':
            raise Exception("Given types are not yet supported... ")
        else:
            raise Exception(
                "Wrong correlation 'type' given. Must be one of 'spec-gcf', 'fam-gcf', ..."
            )

        logger.debug(f"Calculating correlation matrices of type: {type}")

        # Calculate correlation matrix from co-occurence matrices
        M_type1_gcf, M_type1_notgcf, M_nottype1_gcf, M_nottype1_notgcf = calc_correlation_matrix(
            M_type1_strain, self.gcf_strain_occurrence)

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
