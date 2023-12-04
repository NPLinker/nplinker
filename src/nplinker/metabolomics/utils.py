from nplinker.logconfig import LogConfig
from nplinker.strain_collection import StrainCollection
from .molecular_family import MolecularFamily
from .spectrum import Spectrum


logger = LogConfig.getLogger(__name__)


def add_annotation_to_spectrum(annotations: dict[str, dict], spectra: list[Spectrum]) -> None:
    """Add GNPS annotations to the `Spectrum.gnps_annotaions` attribute for input spectra.

    It is possible that some spectra don't have annotations.
    Note that the input `spectra` list is changed in place.

    Args:
        annotations(dict[str, dict]): A dictionary of GNPS annotations, where the keys are
            spectrum ids and the values are GNPS annotations.
        spectra(list[Spectrum]): A list of Spectrum objects.
    """
    for spec in spectra:
        if spec.spectrum_id in annotations:
            spec.gnps_annotations = annotations[spec.spectrum_id]


def add_strains_to_spectrum(
    strains: StrainCollection, spectra: list[Spectrum]
) -> tuple[list[Spectrum], list[Spectrum]]:
    """Add `Strain` objects to the `Spectrum.strains` attribute for input spectra.

    Note that the input `spectra` list is changed in place.

    Args:
        strains(StrainCollection): A collection of strain objects.
        spectra(list[Spectrum]): A list of Spectrum objects.

    Returns:
        tuple(list[Spectrum], list[Spectrum]): A tuple of two lists of Spectrum
            objects. The first list contains Spectrum objects that are updated
            with Strain objects; the second list contains Spectrum objects that
            are not updated with Strain objects becuase no Strain objects are found.
    """
    spectra_with_strains = []
    spectra_without_strains = []
    for spec in spectra:
        try:
            strain_list = strains.lookup(spec.spectrum_id)
        except ValueError:
            spectra_without_strains.append(spec)
            continue

        for strain in strain_list:
            spec.strains.add(strain)
        spectra_with_strains.append(spec)

    logger.info(
        f"{len(spectra_with_strains)} Spectrum objects updated with Strain objects.\n"
        f"{len(spectra_without_strains)} Spectrum objects not updated with Strain objects."
    )

    return spectra_with_strains, spectra_without_strains


def add_spectrum_to_mf(
    spectra: list[Spectrum], mfs: list[MolecularFamily]
) -> tuple[list[MolecularFamily], list[MolecularFamily], dict[MolecularFamily, set[str]]]:
    """Add Spectrum objects to MolecularFamily objects.

    The attribute of `spectra_ids` of MolecularFamily object contains the ids of Spectrum objects.
    These ids are used to find Spectrum objects from the input `spectra` list. The found Spectrum
    objects are added to the `spectra` attribute of MolecularFamily object. It is possible that
    some spectrum ids are not found in the input `spectra` list, and so their Spectrum objects are
    missing in the MolecularFamily object.

    Note that the input `mfs` list is changed in place.

    Args:
        spectra(list[Spectrum]): A list of Spectrum objects.
        mfs(list[MolecularFamily]): A list of MolecularFamily objects.

    Returns:
        tuple(list[MolecularFamily], list[MolecularFamily], dict[MolecularFamily, set[str]]):
            The first list contains MolecularFamily objects that are updated with Spectrum objects.
            The second list contains MolecularFamily objects that are not updated with Spectrum
            objects (all Spectrum objects are missing).
            The dictionary contains MolecularFamily objects as keys and a set of ids of missing
            Spectrum objects as values.
    """
    spec_dict = {spec.spectrum_id: spec for spec in spectra}
    mf_with_spec = []
    mf_without_spec = []
    mf_missing_spec: dict[MolecularFamily, set[str]] = {}
    for mf in mfs:
        for spec_id in mf.spectra_ids:
            try:
                spec = spec_dict[spec_id]
            except KeyError:
                if mf not in mf_missing_spec:
                    mf_missing_spec[mf] = {spec_id}
                else:
                    mf_missing_spec[mf].add(spec_id)
                continue
            mf.add_spectrum(spec)

        if mf.spectra:
            mf_with_spec.append(mf)
        else:
            mf_without_spec.append(mf)

    logger.info(
        f"{len(mf_with_spec)} MolecularFamily objects updated with Spectrum objects.\n"
        f"{len(mf_without_spec)} MolecularFamily objects not updated with Spectrum objects.\n"
        f"{len(mf_missing_spec)} MolecularFamily objects have missing Spectrum objects."
    )
    return mf_with_spec, mf_without_spec, mf_missing_spec


def get_spectra_from_mfs(mfs: list[MolecularFamily]) -> list[Spectrum]:
    """Get all Spectrum objects from given MolecularFamily objects.

    Args:
        mfs(list[MolecularFamily]): A list of MolecularFamily objects.

    Returns:
        list[Spectrum]: A list of Spectrum objects.
    """
    s = set()
    for mf in mfs:
        s |= set(mf.spectra)
    return list(s)
