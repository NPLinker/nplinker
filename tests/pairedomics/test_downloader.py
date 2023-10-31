import os
from pathlib import Path
import httpx
import pytest
from pytest_lazyfixture import lazy_fixture
from nplinker.pairedomics.downloader import PODPDownloader


@pytest.mark.parametrize("expected", [
    Path(os.getenv('HOME'), 'nplinker_data', 'pairedomics'),
    lazy_fixture('tmp_path')
])
def test_default(expected: Path):
    gnps_id = "MSV000079284"

    sut = PODPDownloader(gnps_id, root_dir=str(expected))

    assert sut.gnps_massive_id == gnps_id
    assert sut.working_dir == str(expected)

    assert sut.downloads_dir == str(expected / 'downloads')
    assert sut.project_downloads_dir == str(expected / 'downloads' / gnps_id)

    assert sut.results_dir == str(expected / 'extracted')
    assert sut.project_results_dir == str(expected / 'extracted' / gnps_id)
    assert os.path.exists(str(expected / 'extracted' / gnps_id / 'antismash'))
    assert os.path.exists(str(expected / 'extracted' / gnps_id / 'bigscape'))

    assert sut.all_projects_json_file == str(expected / 'all_projects.json')
    assert sut.project_json_file == str(expected / f"{gnps_id}.json")


def test_download_metabolomics_zipfile(tmp_path):
    sut = PODPDownloader("MSV000079284", root_dir=tmp_path)
    try:
        sut._download_metabolomics_zipfile("c22f44b14a3d450eb836d607cb9521bb")
        expected_path = os.path.join(sut.project_downloads_dir,
                                     'METABOLOMICS-SNETS-c22f44b14a3d450eb836d607cb9521bb.zip')

        assert os.path.exists(expected_path)
        assert (Path(sut.project_results_dir) /
                "molecular_families.tsv").is_file()
        assert (Path(sut.project_results_dir) / "file_mappings.tsv").is_file()
        assert (Path(sut.project_results_dir) / "spectra.mgf").is_file()
    except httpx.TimeoutException:
        pytest.skip("GNPS is down")
