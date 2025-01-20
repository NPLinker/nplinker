from pathlib import Path
import pytest
from nplinker.metabolomics.gnps import GNPSExtractor
from nplinker.metabolomics.gnps import GNPSFormat


#
# Tests for GNPS1 data
#


def test_unknown_workflow_gnps1(gnps_zip_files, tmpdir):
    with pytest.raises(ValueError, match="Unknown workflow type for GNPS archive .*"):
        GNPSExtractor(gnps_zip_files[GNPSFormat.Unknown], tmpdir)


@pytest.mark.parametrize("workflow", [GNPSFormat.FBMN, GNPSFormat.SNETS, GNPSFormat.SNETSV2])
def test_extract_gnps1(
    workflow,
    gnps_zip_files,
    tmpdir,
):
    extractor = GNPSExtractor(gnps_zip_files[workflow], tmpdir)
    assert extractor.gnps_format == workflow
    assert extractor.extract_dir == tmpdir
    assert len(list(Path(tmpdir).iterdir())) == 4

    file_mappings_file = (
        tmpdir / "file_mappings.csv"
        if workflow == GNPSFormat.FBMN
        else tmpdir / "file_mappings.tsv"
    )
    spec_file = tmpdir / "spectra.mgf"
    mf_file = tmpdir / "molecular_families.tsv"
    annotations_file = tmpdir / "annotations.tsv"

    assert file_mappings_file.exists()
    assert spec_file.exists()
    assert mf_file.exists()
    assert annotations_file.exists()


#
# Tests for GNPS2 data
#


def test_unknown_workflow_gnps2(gnps2_tar_files, tmpdir):
    with pytest.raises(ValueError, match="Unknown workflow type for GNPS archive .*"):
        GNPSExtractor(gnps2_tar_files[GNPSFormat.Unknown], tmpdir)


@pytest.mark.parametrize("workflow", [GNPSFormat.GNPS2CN, GNPSFormat.GNPS2FBMN])
def test_extract_gnps2(workflow, gnps2_tar_files, tmpdir):
    extractor = GNPSExtractor(gnps2_tar_files[workflow], tmpdir)
    assert extractor.gnps_format == workflow
    assert extractor.extract_dir == tmpdir
    assert len(list(Path(tmpdir).iterdir())) == 4

    file_mappings_file = tmpdir / "file_mappings.csv"
    spec_file = tmpdir / "spectra.mgf"
    mf_file = tmpdir / "molecular_families.tsv"
    annotations_file = tmpdir / "annotations.tsv"

    assert file_mappings_file.exists()
    assert spec_file.exists()
    assert mf_file.exists()
    assert annotations_file.exists()
