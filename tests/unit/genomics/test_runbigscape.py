import pytest
from nplinker.genomics import bigscape
from .. import DATA_DIR


class TestRunBigscape:
    def test_run_bigscape_v1(self, tmp_path):
        result = bigscape.run_bigscape(
            antismash_path=DATA_DIR,
            output_path=tmp_path,
            extra_params="--help",
            version=1,
        )

        assert result is True

    def test_run_bigscape_v2(self, tmp_path):
        result = bigscape.run_bigscape(
            antismash_path=DATA_DIR,
            output_path=tmp_path,
            extra_params="--help",
            version=2,
        )

        assert result is True

    def test_run_bigscape_small_dataset_v1(self, tmp_path):
        result = bigscape.run_bigscape(
            antismash_path=DATA_DIR / "bigscape/minimal_dataset",
            output_path=tmp_path,
            extra_params="",
            version=1,
        )

        assert result is True

    def test_run_bigscape_small_dataset_v2(self, tmp_path):
        result = bigscape.run_bigscape(
            antismash_path=DATA_DIR / "bigscape/minimal_dataset",
            output_path=tmp_path,
            extra_params="",
            version=2,
        )

        assert result is True

    def test_run_bigscape_wrong_version(self, tmp_path):
        with pytest.raises(ValueError) as e:
            bigscape.run_bigscape(
                antismash_path=DATA_DIR,
                output_path=tmp_path,
                extra_params="--help",
                version=3,
            )

        assert "version" in e.value.args[0]

    def test_input_path_not_exist_v1(self, tmp_path):
        with pytest.raises(FileNotFoundError) as e:
            bigscape.run_bigscape(
                antismash_path=tmp_path / "not_exist",
                output_path=tmp_path,
                extra_params="",
                version=1,
            )

        assert "antismash_path" in e.value.args[0]

    def test_input_path_not_exist_v2(self, tmp_path):
        with pytest.raises(FileNotFoundError) as e:
            bigscape.run_bigscape(
                antismash_path=tmp_path / "not_exist",
                output_path=tmp_path,
                extra_params="",
                version=2,
            )

        assert "antismash_path" in e.value.args[0]
