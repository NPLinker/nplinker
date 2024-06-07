from __future__ import annotations
from typing import TYPE_CHECKING
from typing import Sequence
import numpy as np
import pandas as pd


if TYPE_CHECKING:
    from nplinker.genomics import GCF
    from nplinker.metabolomics import MolecularFamily
    from nplinker.metabolomics import Spectrum
    from nplinker.strain import StrainCollection


def isinstance_all(*objects, objtype) -> bool:
    """Check if all objects are of the given type."""
    return all(isinstance(x, objtype) for x in objects)


def get_presence_gcf_strain(gcfs: Sequence[GCF], strains: StrainCollection) -> pd.DataFrame:
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


def get_presence_spec_strain(
    spectra: Sequence[Spectrum], strains: StrainCollection
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


def get_presence_mf_strain(
    mfs: Sequence[MolecularFamily], strains: StrainCollection
) -> pd.DataFrame:
    """Get the occurence of strains in molecular families.

    The occurence is a DataFrame with molecular families as rows and
    strains as columns, where index is `mf.id` and column name is
    `strain.id`. The values are 1 if the molecular family contains the
    strain and 0 otherwise.
    """
    df_mf_strain = pd.DataFrame(
        np.zeros((len(mfs), len(strains))),
        index=[mf.id for mf in mfs],
        columns=[strain.id for strain in strains],
        dtype=int,
    )
    for mf in mfs:
        for strain in strains:
            if mf.has_strain(strain):
                df_mf_strain.loc[mf.id, strain.id] = 1
    return df_mf_strain
