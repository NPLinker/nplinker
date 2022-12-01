from pathlib import Path
from shutil import rmtree
from tempfile import mkdtemp
import zipfile
import pytest
from nplinker import utils
from nplinker.utils import find_delimiter
from tests import DATA_DIR


@pytest.mark.parametrize(
    "filename, expected",
    [[DATA_DIR / "nodes.tsv", '\t'], [DATA_DIR / "nodes_mwe.csv", ',']])
def test_find_delimiter(filename, expected):
    actual = find_delimiter(filename)
    assert actual == expected


BGC_GBK_URL = "https://mibig.secondarymetabolites.org/repository/BGC0000001/BGC0000001.gbk"
MIBIG_METADATAS_URL = "https://dl.secondarymetabolites.org/mibig/mibig_json_3.1.tar.gz"


class TestDownloadUrl:
    """Test utils.download_url"""

    @pytest.fixture
    def temppath1(self):
        temppath = mkdtemp()
        yield temppath
        rmtree(temppath)

    def test_default(self, temppath1):
        utils.download_url(url=BGC_GBK_URL, root=temppath1)
        f = Path(temppath1) / "BGC0000001.gbk"
        assert f.is_file()

    def test_optional_args(self, temppath1):
        utils.download_url(url=BGC_GBK_URL,
                           root=temppath1,
                           filename="example.gbk",
                           md5="f33f98481e6d7855c7622d715ef148ac")
        f = Path(temppath1) / "example.gbk"
        assert f.is_file()


class TestExtractArchive:
    """Test utils.extract_archive"""

    @pytest.fixture
    def archive(self):
        temppath = mkdtemp()
        utils.download_url(MIBIG_METADATAS_URL, temppath)
        archive = Path(temppath) / "mibig_json_3.1.tar.gz"
        yield archive

    def test_default(self, archive):
        utils.extract_archive(archive)
        dir = archive.parent / "mibig_json_3.1"
        assert dir.exists()
        assert dir.is_dir()

    def test_optional_args(self, archive):
        to_path = mkdtemp()
        utils.extract_archive(archive, to_path=to_path, remove_finished=True)
        dir = Path(to_path) / "mibig_json_3.1"
        assert dir.exists()
        assert dir.is_dir()
        assert not archive.exists()

class TestDownloadAndExtractArchive:
    """Test utils.download_and_extract_archive"""

    @pytest.fixture
    def temppath1(self):
        temppath = mkdtemp()
        yield temppath
        rmtree(temppath)

    @pytest.fixture
    def temppath2(self):
        temppath = mkdtemp()
        yield temppath
        rmtree(temppath)

    def test_defaults(self, temppath1):
        utils.download_and_extract_archive(url=MIBIG_METADATAS_URL,
                                           download_root=temppath1)

        fdownload = Path(temppath1) / "mibig_json_3.1.tar.gz"
        fextract = Path(temppath1) / "mibig_json_3.1"

        assert fdownload.is_file()
        assert fextract.is_dir()

    def test_optional_args(self, temppath1, temppath2):
        utils.download_and_extract_archive(
            url=MIBIG_METADATAS_URL,
            download_root=temppath1,
            extract_root=temppath2,
            filename="example.tar.gz",
            md5="643d1349722a9437d8dcf558dac5f815")

        fdownload = Path(temppath1) / "example.tar.gz"
        fextract = Path(temppath2) / "mibig_json_3.1"

        assert fdownload.is_file()
        assert fextract.exists()
        assert fextract.is_dir()

    def test_arg_remove_finished(self, temppath1):
        utils.download_and_extract_archive(url=MIBIG_METADATAS_URL,
                                           download_root=temppath1,
                                           remove_finished=True)

        fdownload = Path(temppath1) / "mibig_json_3.1.tar.gz"
        fextract = Path(temppath1) / "mibig_json_3.1"

        assert not fdownload.is_file()
        assert fextract.is_dir()


class TestListDirs:
    """Test utils.list_dirs"""

    root = Path(__file__).parent

    def test_default(self):
        dirs = utils.list_dirs(root=self.root)
        assert isinstance(dirs, list)
        assert len(dirs) >= 1
        assert "data" in dirs
        assert str(self.root / "data") not in dirs

    def test_prefix(self):
        dirs = utils.list_dirs(root=self.root, prefix=True)
        assert isinstance(dirs, list)
        assert len(dirs) >= 1
        assert "data" not in dirs
        assert str(self.root / "data") in dirs


class TestListFiles:
    """Test utils.list_files"""

    root = Path(__file__).parent

    def test_default(self):
        files = utils.list_files(self.root, suffix=".py")
        assert isinstance(files, list)
        assert "test_utils.py" in files
        assert str(self.root / "test_utils.py") not in files

    def test_prefix(self):
        files = utils.list_files(self.root, suffix=".py", prefix=True)
        assert isinstance(files, list)
        assert "test_utils.py" not in files
        assert str(self.root / "test_utils.py") in files


def test_extract_file_matching_patter(tmp_path):
    archive = zipfile.ZipFile(DATA_DIR / "ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-92036537-download_cytoscape_data.zip")
    utils.extract_file_matching_pattern(
        archive,
        "quantification_table_reformatted",
        ".csv" ,
        tmp_path,
        "test.csv"
    )
    expected : Path = tmp_path / "test.csv"
    
    assert expected.exists()
    assert expected.is_file()