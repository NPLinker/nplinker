import filecmp
from pathlib import Path
import pytest
from nplinker.metabolomics.gnps.gnps_extractor import GNPSExtractor
from nplinker.metabolomics.gnps.gnps_format import GNPSFormat


def test_unknown_workflow(gnps_zip_files, tmpdir):
    with pytest.raises(ValueError) as e:
        GNPSExtractor(gnps_zip_files[GNPSFormat.Unknown], tmpdir)
    assert "Unknown workflow type for GNPS archive" in str(e.value)


@pytest.mark.parametrize(
    "workflow", [GNPSFormat.FBMN, GNPSFormat.SNETS, GNPSFormat.SNETSV2])
def test_supported_workflows(workflow, gnps_zip_files, tmpdir):
    extractor = GNPSExtractor(gnps_zip_files[workflow], tmpdir)
    assert extractor.gnps_format == workflow
    assert extractor.extract_dir == tmpdir


@pytest.mark.parametrize(
    "workflow", [GNPSFormat.FBMN, GNPSFormat.SNETS, GNPSFormat.SNETSV2])
def test_extract(workflow, gnps_zip_files, gnps_file_mappings_files,
                 gnps_spectra_files, gnps_mf_files, gnps_annotations_files,
                 tmpdir):
    GNPSExtractor(gnps_zip_files[workflow], tmpdir)
    assert len(list(Path(tmpdir).iterdir())) == 4

    file_mappings_file = tmpdir / "file_mappings.csv" if workflow == GNPSFormat.FBMN else tmpdir / "file_mappings.tsv"
    spec_file = tmpdir / "spectra.mgf"
    mf_file = tmpdir / "molecular_families.tsv"
    annotations_file = tmpdir / "annotations.tsv"

    assert file_mappings_file.exists()
    assert spec_file.exists()
    assert mf_file.exists()
    assert annotations_file.exists()

    assert filecmp.cmp(file_mappings_file, gnps_file_mappings_files[workflow])
    assert filecmp.cmp(spec_file, gnps_spectra_files[workflow])
    assert filecmp.cmp(mf_file, gnps_mf_files[workflow])
    assert filecmp.cmp(annotations_file, gnps_annotations_files[workflow])
