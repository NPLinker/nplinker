from __future__ import annotations
from collections.abc import Sequence
import pandas as pd
from nplinker.genomics import GCF
from nplinker.metabolomics import MolecularFamily
from nplinker.metabolomics import Spectrum
from nplinker.strain import StrainCollection


def get_presence_gcf_strain(gcfs: Sequence[GCF], strains: StrainCollection) -> pd.DataFrame:
    """Get the occurrence of strains in gcfs.

    The occurrence is a DataFrame with GCF objects as index and Strain objects as columns, and the
    values are 1 if the gcf occurs in the strain,  0 otherwise.
    """
    df_gcf_strain = pd.DataFrame(
        0,
        index=gcfs,
        columns=list(strains),
        dtype=int,
    )  # type: ignore
    for gcf in gcfs:
        for strain in strains:
            if gcf.has_strain(strain):
                df_gcf_strain.loc[gcf, strain] = 1
    return df_gcf_strain  # type: ignore


def get_presence_spec_strain(
    spectra: Sequence[Spectrum], strains: StrainCollection
) -> pd.DataFrame:
    """Get the occurrence of strains in spectra.

    The occurrence is a DataFrame with Spectrum objects as index and Strain objects as columns, and
    the values are 1 if the spectrum occurs in the strain, 0 otherwise.
    """
    df_spec_strain = pd.DataFrame(
        0,
        index=spectra,
        columns=list(strains),
        dtype=int,
    )  # type: ignore
    for spectrum in spectra:
        for strain in strains:
            if spectrum.has_strain(strain):
                df_spec_strain.loc[spectrum, strain] = 1
    return df_spec_strain  # type: ignore


def get_presence_mf_strain(
    mfs: Sequence[MolecularFamily], strains: StrainCollection
) -> pd.DataFrame:
    """Get the occurrence of strains in molecular families.

    The occurrence is a DataFrame with MolecularFamily objects as index and Strain objects as
    columns, and the values are 1 if the molecular family occurs in the strain, 0 otherwise.
    """
    df_mf_strain = pd.DataFrame(
        0,
        index=mfs,
        columns=list(strains),
        dtype=int,
    )  # type: ignore
    for mf in mfs:
        for strain in strains:
            if mf.has_strain(strain):
                df_mf_strain.loc[mf, strain] = 1
    return df_mf_strain  # type: ignore
