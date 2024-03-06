from __future__ import annotations
from typing import TYPE_CHECKING
from typing import Sequence
import numpy as np
import pandas as pd
from nplinker.genomics.gcf import GCF
from nplinker.logconfig import LogConfig
from nplinker.metabolomics import MolecularFamily
from nplinker.metabolomics import Spectrum
from .utils import calc_correlation_matrix
from .utils import isinstance_all


if TYPE_CHECKING:
    from nplinker.strain import Strain
    from nplinker.strain import StrainCollection

logger = LogConfig.getLogger(__name__)

LINK_TYPES = ["spec-gcf", "mf-gcf"]


class DataLinks:
    def __init__(
        self,
        gcfs: Sequence[GCF],
        spectra: Sequence[Spectrum],
        mfs: Sequence[MolecularFamily],
        strains: StrainCollection,
    ):
        """DataLinks class to store occurrence and co-occurrence information.

        Occurrence refers to the presence of a spectrum/gcf/mf in a strain.
        Co-occurrence refers to the presence of a spectrum/mf and a gcf in a strain.

        Args:
            gcfs: A list of GCF objects.
            spectra: A list of Spectrum objects.
            mfs: A list of MolecularFamily objects.
            strains: A StrainCollection object.

        Attributes:
            occurrence_gcf_strain: A DataFrame to store occurrence of
                gcfs with respect to strains.
            occurrence_spec_strain: A DataFrame to store occurrence of
                spectra with respect to strains.
            occurrence_mf_strain: A DataFrame to store occurrence of
                molecular families with respect to strains.
            cooccurrence_spec_gcf: A DataFrame to store co-occurrence
                of the presence of spectra and the presence of gcfs with respect
                to strains.
            cooccurrence_spec_notgcf: A DataFrame to store co-occurrence
                of the presence of spectra and the absence of gcfs with respect
                to strains. "notgcf" means the absence of gcfs.
            cooccurrence_notspec_gcf: A DataFrame to store co-occurrence
                of the absence of spectra and the presence of gcfs with respect
                to strains. "notspec" means the absence of spectra.
            cooccurrence_notspec_notgcf: A DataFrame to store co-occurrence
                of the absence of spectra and the absence of gcfs with respect
                to strains.
            cooccurrence_mf_gcf: A DataFrame to store co-occurrence
                of the presence of molecular families and the presence of gcfs
                with respect to strains.
            cooccurrence_mf_notgcf: A DataFrame to store co-occurrence
                of the presence of molecular families and the absence of gcfs
                with respect to strains. "notgcf" means the absence of gcfs.
            cooccurrence_notmf_gcf: A DataFrame to store co-occurrence
                of the absence of molecular families and the presence of gcfs
                with respect to strains. "notmf" means the absence of molecular
                families.
            cooccurrence_notmf_notgcf: A DataFrame to store co-occurrence
                of the absence of molecular families and the absence of gcfs
                with respect to strains.
        """
        self._strains = strains

        logger.debug(
            "Create occurrence dataframes: spectra<->strains, gcfs<->strains and mfs<->strains."
        )
        # DataFrame to store occurrence of gcfs/spectra/mfs with respect to strains
        # values = 1 where gcf/spec/fam occur in strain, 0 otherwise
        self.occurrence_gcf_strain = self._get_occurrence_gcf_strain(gcfs, strains)
        self.occurrence_spec_strain = self._get_occurrence_spec_strain(spectra, strains)
        self.occurrence_mf_strain = self._get_occurrence_mf_strain(mfs, strains)

        # DataFrame to store co-occurrence of "spectra<->gcf" or "mfs<->gcf"
        logger.debug("Create correlation matrices: spectra<->gcfs.")
        (
            self.cooccurrence_spec_gcf,
            self.cooccurrence_spec_notgcf,
            self.cooccurrence_notspec_gcf,
            self.cooccurrence_notspec_notgcf,
        ) = self._get_cooccurrence(link_type="spec-gcf")
        logger.debug("Create correlation matrices: mol-families<->gcfs.")
        (
            self.cooccurrence_mf_gcf,
            self.cooccurrence_mf_notgcf,
            self.cooccurrence_notmf_gcf,
            self.cooccurrence_notmf_notgcf,
        ) = self._get_cooccurrence(link_type="mf-gcf")

    def get_common_strains(
        self,
        spectra_or_mfs: Sequence[Spectrum | MolecularFamily],
        gcfs: Sequence[GCF],
        filter_no_shared: bool = False,
    ) -> dict[tuple[Spectrum | MolecularFamily, GCF], list[Strain]]:
        """Get common strains between given spectra/molecular families and GCFs.

        Args:
            spectra_or_mfs:
                A list of Spectrum and/or MolecularFamily objects.
            gcfs: A list of GCF objects.
            filter_no_shared: If True, the pairs of spectrum/mf and GCF
                without common strains will be removed from the returned dict;

        Returns:
            Keys are tuples of (Spectrum/MolecularFamily, GCF) and values are a list of shared
            Strain objects.

        Raises:
            ValueError: If given `spectra_or_mfs` or `gcfs` is empty.
            TypeError: If given `spectra_or_mfs` or `gcfs` is not a list of
                Spectrum/MolecularFamily or GCF objects, respectively.
        """
        # Check input arguments
        if len(spectra_or_mfs) == 0 or len(gcfs) == 0:
            raise ValueError("Empty list for first or second argument.")
        if not isinstance_all(*spectra_or_mfs, objtype=(Spectrum, MolecularFamily)):
            raise TypeError("First argument must be Spectrum and/or MolecularFamily objects.")
        if not isinstance_all(*gcfs, objtype=GCF):
            raise TypeError("Second argument must be GCF objects.")

        # Assume that 3 occurrence dataframes have same df.columns (strain ids)
        strain_ids = self.occurrence_gcf_strain.columns
        results = {}
        for obj in spectra_or_mfs:
            for gcf in gcfs:
                if isinstance(obj, Spectrum):
                    shared_strains = strain_ids[
                        np.logical_and(
                            self.occurrence_spec_strain.loc[obj.spectrum_id],
                            self.occurrence_gcf_strain.loc[gcf.gcf_id],
                        )
                    ]
                else:
                    shared_strains = strain_ids[
                        np.logical_and(
                            self.occurrence_mf_strain.loc[obj.family_id],
                            self.occurrence_gcf_strain.loc[gcf.gcf_id],
                        )
                    ]
                if filter_no_shared and len(shared_strains) == 0:
                    continue
                results[(obj, gcf)] = [
                    strain
                    for strain_id in shared_strains
                    for strain in self._strains.lookup(strain_id)
                ]
        return results

    def _get_occurrence_gcf_strain(
        self, gcfs: Sequence[GCF], strains: StrainCollection
    ) -> pd.DataFrame:
        """Get the occurence of strains in gcfs.

        The occurence is a DataFrame with gcfs as rows and strains as columns,
        where index is `gcf.gcf_id` and column name is `strain.id`. The values
        are 1 if the gcf contains the strain and 0 otherwise.
        """
        df_gcf_strain = pd.DataFrame(
            np.zeros((len(gcfs), len(strains))),
            index=[gcf.gcf_id for gcf in gcfs],
            columns=[strain.id for strain in strains],
            dtype=int,
        )
        for gcf in gcfs:
            for strain in strains:
                if gcf.has_strain(strain):
                    df_gcf_strain.loc[gcf.gcf_id, strain.id] = 1
        return df_gcf_strain

    def _get_occurrence_spec_strain(
        self, spectra: Sequence[Spectrum], strains: StrainCollection
    ) -> pd.DataFrame:
        """Get the occurence of strains in spectra.

        The occurence is a DataFrame with spectra as rows and strains as columns,
        where index is `spectrum.spectrum_id` and column name is `strain.id`.
        The values are 1 if the spectrum contains the strain and 0 otherwise.
        """
        df_spec_strain = pd.DataFrame(
            np.zeros((len(spectra), len(strains))),
            index=[spectrum.spectrum_id for spectrum in spectra],
            columns=[strain.id for strain in strains],
            dtype=int,
        )
        for spectrum in spectra:
            for strain in strains:
                if spectrum.has_strain(strain):
                    df_spec_strain.loc[spectrum.spectrum_id, strain.id] = 1
        return df_spec_strain

    def _get_occurrence_mf_strain(
        self, mfs: Sequence[MolecularFamily], strains: StrainCollection
    ) -> pd.DataFrame:
        """Get the occurence of strains in molecular families.

        The occurence is a DataFrame with molecular families as rows and
        strains as columns, where index is `mf.family_id` and column name is
        `strain.id`. The values are 1 if the molecular family contains the
        strain and 0 otherwise.
        """
        df_mf_strain = pd.DataFrame(
            np.zeros((len(mfs), len(strains))),
            index=[mf.family_id for mf in mfs],
            columns=[strain.id for strain in strains],
            dtype=int,
        )
        for mf in mfs:
            for strain in strains:
                if mf.has_strain(strain):
                    df_mf_strain.loc[mf.family_id, strain.id] = 1
        return df_mf_strain

    def _get_cooccurrence(
        self, link_type: str = "spec-gcf"
    ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """Calculate co-occurrence for given link types across strains.

        Args:
            link_type: Type of link to calculate co-occurrence for,
                either 'spec-gcf' or 'mf-gcf'.
        """
        if link_type == "spec-gcf":
            met_strain_occurrence = self.occurrence_spec_strain
        elif link_type == "mf-gcf":
            met_strain_occurrence = self.occurrence_mf_strain
        else:
            raise ValueError(
                f"Link type {link_type} is not supported. Use 'spec-gcf' or 'mf-gcf' instead."
            )
        logger.debug(f"Calculating correlation matrices of type: {link_type}")
        m1, m2, m3, m4 = calc_correlation_matrix(met_strain_occurrence, self.occurrence_gcf_strain)
        df_met_gcf = pd.DataFrame(
            m1,
            index=met_strain_occurrence.index,
            columns=self.occurrence_gcf_strain.index,
            dtype=int,
        )
        df_met_notgcf = pd.DataFrame(
            m2,
            index=met_strain_occurrence.index,
            columns=self.occurrence_gcf_strain.index,
            dtype=int,
        )
        df_notmet_gcf = pd.DataFrame(
            m3,
            index=met_strain_occurrence.index,
            columns=self.occurrence_gcf_strain.index,
            dtype=int,
        )
        df_notmet_notgcf = pd.DataFrame(
            m4,
            index=met_strain_occurrence.index,
            columns=self.occurrence_gcf_strain.index,
            dtype=int,
        )
        return df_met_gcf, df_met_notgcf, df_notmet_gcf, df_notmet_notgcf
