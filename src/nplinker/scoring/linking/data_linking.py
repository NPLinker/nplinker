from __future__ import annotations
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
    from nplinker.strains import Strain

logger = LogConfig.getLogger(__name__)

LINK_TYPES = ['spec-gcf', 'mf-gcf']


class DataLinks():

    def __init__(self, gcfs: Sequence[GCF], spectra: Sequence[Spectrum],
                 mfs: Sequence[MolecularFamily], strains: StrainCollection):
        """DataLinks class to store occurrence and co-occurrence information.

        Occurrence refers to the presence of a spectrum/gcf/mf in a strain.
        Co-occurrence refers to the presence of a spectrum/mf and a gcf in a strain.

        Args:
            gcfs(Sequence[GCF]): A list of GCF objects.
            spectra(Sequence[Spectrum]): A list of Spectrum objects.
            mfs(Sequence[MolecularFamily]): A list of MolecularFamily objects.
            strains(StrainCollection): A StrainCollection object.

        Attributes:
            occurrence_gcf_strain(pd.DataFrame): A DataFrame to store occurrence of
                gcfs with respect to strains.
            occurrence_spec_strain(pd.DataFrame): A DataFrame to store occurrence of
                spectra with respect to strains.
            occurrence_mf_strain(pd.DataFrame): A DataFrame to store occurrence of
                molecular families with respect to strains.
            cooccurrence_spec_gcf(pd.DataFrame): A DataFrame to store co-occurrence
                of spectra<->gcfs.
            cooccurrence_spec_notgcf(pd.DataFrame): A DataFrame to store co-occurrence
                of spectra<->not gcfs.
            cooccurrence_notspec_gcf(pd.DataFrame): A DataFrame to store co-occurrence
                of not spectra<->gcfs.
            cooccurrence_notspec_notgcf(pd.DataFrame): A DataFrame to store co-occurrence
                of not spectra<->not gcfs.
            cooccurrence_mf_gcf(pd.DataFrame): A DataFrame to store co-occurrence
                of molecular families<->gcfs.
            cooccurrence_mf_notgcf(pd.DataFrame): A DataFrame to store co-occurrence
                of molecular families<->not gcfs.
            cooccurrence_notmf_gcf(pd.DataFrame): A DataFrame to store co-occurrence
                of not molecular families<->gcfs.
            cooccurrence_notmf_notgcf(pd.DataFrame): A DataFrame to store co-occurrence
                of not molecular families<->not gcfs.
            mapping_gcf(pd.DataFrame): A DataFrame to store mappings for gcfs.
            mapping_spec(pd.DataFrame): A DataFrame to store mappings for spectra.
            mapping_mf(pd.DataFrame): A DataFrame to store mappings for molecular families.
            mapping_strain(pd.DataFrame): A DataFrame to store mappings for strains.
        """
        self._strains = strains

        logger.debug(
            "Create occurrence dataframes: spectra<->strains, gcfs<->strains and mfs<->strains."
        )
        # DataFrame to store occurrence of gcfs/spectra/mfs with respect to strains
        # values = 1 where gcf/spec/fam occur in strain, 0 otherwise
        self.occurrence_gcf_strain = self._get_occurrence_gcf_strain(
            gcfs, strains)
        self.occurrence_spec_strain = self._get_occurrence_spec_strain(
            spectra, strains)
        self.occurrence_mf_strain = self._get_occurrence_mf_strain(
            mfs, strains)

        # DataFrame to store mapping tables, check `_get_mappings_from_occurance` for details
        # TODO: these mappings could be removed when refactoring LinkFinder
        self.mapping_spec = pd.DataFrame()
        self.mapping_gcf = pd.DataFrame()
        self.mapping_fam = pd.DataFrame()
        self.mapping_strain = pd.DataFrame()
        self._get_mappings_from_occurrence()

        # DataFrame to store co-occurrence of "spectra<->gcf" or "mfs<->gcf"
        logger.debug("Create correlation matrices: spectra<->gcfs.")
        (self.cooccurrence_spec_gcf, self.cooccurrence_spec_notgcf,
         self.cooccurrence_notspec_gcf,
         self.cooccurrence_notspec_notgcf) = self._get_cooccurrence(
             link_type='spec-gcf')
        logger.debug("Create correlation matrices: mol-families<->gcfs.")
        (self.cooccurrence_mf_gcf, self.cooccurrence_mf_notgcf,
         self.cooccurrence_notmf_gcf,
         self.cooccurrence_notmf_notgcf) = self._get_cooccurrence(
             link_type='mf-gcf')

    def get_common_strains(
        self,
        spectra_or_mfs: Sequence[Spectrum] | Sequence[MolecularFamily],
        gcfs: Sequence[GCF],
        filter_no_shared: bool = False
    ) -> dict[tuple[Spectrum | MolecularFamily, GCF], list[Strain]]:
        """Get common strains between given spectra/molecular families and GCFs.

        Note that SingletonFamily objects are excluded from given `spectra_or_mfs`.

        Args:
            spectra_or_mfs(Sequence[Spectrum] | Sequence[MolecularFamily]):
                A list of Spectrum or MolecularFamily objects.
            gcfs(Sequence[GCF]): A list of GCF objects.
            filter_no_shared(bool): If True, the pairs of spectrum/mf and GCF
                without common strains will be removed from the returned dict;

        Returns:
            dict: A dict where the keys are tuples of (Spectrum/MolecularFamily, GCF)
            and values are a list of shared Strain objects.
        """
        if not isinstance(spectra_or_mfs[0], (Spectrum, MolecularFamily)):
            raise ValueError(
                'Must provide Spectra or MolecularFamilies as the first argument!'
            )
        if not isinstance(gcfs[0], GCF):
            raise ValueError('Must provide GCFs as the second argument!')

        # Assume that 3 occurrence dataframes have same df.columns (strain ids)
        strain_ids = self.occurrence_gcf_strain.columns
        results = {}
        for obj in spectra_or_mfs:
            if isinstance(obj, SingletonFamily):
                continue
            for gcf in gcfs:
                if isinstance(obj, Spectrum):
                    shared_strains = strain_ids[np.logical_and(
                        self.occurrence_spec_strain.loc[obj.spectrum_id],
                        self.occurrence_gcf_strain.loc[gcf.gcf_id])]
                else:
                    shared_strains = strain_ids[np.logical_and(
                        self.occurrence_mf_strain.loc[obj.family_id],
                        self.occurrence_gcf_strain.loc[gcf.gcf_id])]
                if filter_no_shared and len(shared_strains) == 0:
                    continue
                results[(obj, gcf)] = [
                    self._strains.lookup(strain_id)
                    for strain_id in shared_strains
                ]
        return results

    def _get_occurrence_gcf_strain(self, gcfs: Sequence[GCF],
                                   strains: StrainCollection) -> pd.DataFrame:
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
        return df_gcf_strain

    def _get_occurrence_spec_strain(self, spectra: Sequence[Spectrum],
                                    strains: StrainCollection) -> pd.DataFrame:
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
        return df_spec_strain

    def _get_occurrence_mf_strain(self, mfs: Sequence[MolecularFamily],
                                  strains: StrainCollection) -> pd.DataFrame:
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
        return df_mf_strain

    def _get_mappings_from_occurrence(self):

        # pd.Series with index = gcf.gcf_id and value = number of strains where gcf occurs
        self.mapping_gcf["no of strains"] = self.occurrence_gcf_strain.sum(
            axis=1)
        # pd.Series with index = spectrum.spectrum_id and value = number of strains where spec occurs
        self.mapping_spec["no of strains"] = self.occurrence_spec_strain.sum(
            axis=1)
        # pd.Series with index = mf.family_id and value = number of strains where mf occurs
        self.mapping_fam["no of strains"] = self.occurrence_mf_strain.sum(
            axis=1)
        # pd.Series with index = strain.id and value = number of spectra in strain
        self.mapping_strain["no of spectra"] = self.occurrence_spec_strain.sum(
            axis=0)

    def _get_cooccurrence(
        self,
        link_type: str = 'spec-gcf'
    ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """Calculate co-occurrence for given link types across strains.

        Args:
            link_type(str): Type of link to calculate co-occurrence for,
                either 'spec-gcf' or 'mf-gcf'.
        """
        if link_type == 'spec-gcf':
            met_strain_occurrence = self.occurrence_spec_strain
        elif link_type == 'mf-gcf':
            met_strain_occurrence = self.occurrence_mf_strain
        else:
            raise ValueError(
                f"Link type {link_type} is not supported. Use 'spec-gcf' or 'mf-gcf' instead."
            )
        logger.debug(f"Calculating correlation matrices of type: {link_type}")
        m1, m2, m3, m4 = calc_correlation_matrix(met_strain_occurrence,
                                                 self.occurrence_gcf_strain)
        df_met_gcf = pd.DataFrame(m1,
                                  index=met_strain_occurrence.index,
                                  columns=self.occurrence_gcf_strain.index)
        df_met_notgcf = pd.DataFrame(m2,
                                     index=met_strain_occurrence.index,
                                     columns=self.occurrence_gcf_strain.index)
        df_notmet_gcf = pd.DataFrame(m3,
                                     index=met_strain_occurrence.index,
                                     columns=self.occurrence_gcf_strain.index)
        df_notmet_notgcf = pd.DataFrame(
            m4,
            index=met_strain_occurrence.index,
            columns=self.occurrence_gcf_strain.index)
        return df_met_gcf, df_met_notgcf, df_notmet_gcf, df_notmet_notgcf
