import shutil
import tempfile
from pathlib import Path
import pytest
from nplinker.genomics import mibig


class TestDownloadAndExtractMibigMetadatas():

    @pytest.fixture
    def temppath1(self):
        temppath = tempfile.mkdtemp()
        yield temppath
        shutil.rmtree(temppath)

    @pytest.fixture
    def temppath2(self):
        temppath = tempfile.mkdtemp()
        yield temppath
        shutil.rmtree(temppath)

    def test_default(self, temppath1, temppath2):
        mibig.download_and_extract_mibig_metadata(temppath1, temppath2)
        archive = Path(temppath1) / "mibig_json_3.1.tar.gz"
        metadata = Path(temppath2) / "BGC0000002.json"
        assert archive.exists()
        assert archive.is_file()
        assert metadata.exists()
        assert metadata.is_file()

    def test_version(self, temppath1, temppath2):
        mibig.download_and_extract_mibig_metadata(temppath1,
                                                  temppath2,
                                                  version="1.4")
        archive = Path(temppath1) / "mibig_json_1.4.tar.gz"
        metadata = Path(temppath2) / "BGC0000002.json"
        assert archive.exists()
        assert archive.is_file()
        assert metadata.exists()
        assert metadata.is_file()

    def test_error_same_path(self, temppath1):
        with pytest.raises(ValueError) as e:
            mibig.download_and_extract_mibig_metadata(temppath1, temppath1)
        assert e.value.args[
            0] == "Identical path of download directory and extract directory"

    def test_error_nonempty_path(self, temppath1):
        extract_path = Path(__file__).parent
        with pytest.raises(ValueError) as e:
            mibig.download_and_extract_mibig_metadata(temppath1, extract_path)
        assert "Nonempty directory" in e.value.args[0]
