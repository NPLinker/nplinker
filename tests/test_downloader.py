import os
from pathlib import Path

from nplinker.pairedomics.downloader import Downloader

def test_default():
    gnps_id = "MSV000079284"
    sut = Downloader(gnps_id)

    expected = Path(os.getenv('HOME'), 'nplinker_data',
                                        'pairedomics')

    assert sut is not None
    assert sut.gnps_massive_id == gnps_id
    assert sut.local_cache == str(expected)

    assert sut.local_download_cache == str(expected / 'downloads')
    assert sut.project_download_cache == str(expected / 'downloads' / gnps_id)

    assert sut.local_file_cache == str(expected / 'extracted')
    assert sut.project_file_cache == str(expected / 'extracted'/ gnps_id)
    assert sut.strain_mappings_file == str(expected / 'extracted'/ gnps_id / 'strain_mappings.csv')
    assert os.path.exists(str(expected / 'extracted'/ gnps_id / 'antismash'))
    assert os.path.exists(str(expected / 'extracted'/ gnps_id / 'bigscape'))

    assert sut.all_project_json_file == str(expected / 'all_projects.json')
    assert sut.project_json_file == str(expected / f"{gnps_id}.json")


