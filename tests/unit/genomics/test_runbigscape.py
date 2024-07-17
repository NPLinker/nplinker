import os
import pytest
from nplinker.genomics import bigscape
from .. import DATA_DIR


@pytest.mark.parametrize("version", [1, 2])
def test_run_bigscape(tmp_path, version):
    """Test whether BiG-SCAPE runs at all using the --help command"""
    result = bigscape.run_bigscape(
        antismash_path=tmp_path,
        output_path=tmp_path,
        extra_params="--help",
        version=version,
    )

    assert result is True


@pytest.mark.parametrize("version", [1, 2])
def test_run_bigscape_small_dataset(tmp_path, version):
    pytest.skip("This test is too slow to run in CI")
    result = bigscape.run_bigscape(
        antismash_path=DATA_DIR / "bigscape" / "minimal_dataset",
        output_path=tmp_path,
        extra_params="",
        version=version,
    )

    assert result is True


def test_run_bigscape_wrong_version(tmp_path):
    with pytest.raises(ValueError) as e:
        bigscape.run_bigscape(
            antismash_path=DATA_DIR,
            output_path=tmp_path,
            extra_params="--help",
            version=3,
        )

    assert "version" in e.value.args[0]


@pytest.mark.parametrize("version", [1, 2])
def test_input_path_not_exist(tmp_path, version):
    with pytest.raises(FileNotFoundError) as e:
        bigscape.run_bigscape(
            antismash_path=tmp_path / "not_exist",
            output_path=tmp_path,
            extra_params="",
            version=version,
        )

    assert "antismash_path" in e.value.args[0]


@pytest.mark.parametrize("version", [1, 2])
def test_bad_parameters(tmp_path, version):
    with pytest.raises(RuntimeError) as e:
        bigscape.run_bigscape(
            antismash_path=tmp_path,
            output_path=tmp_path,
            extra_params="--this-is-not-a-real-argument",
            version=version,
        )

    assert "BiG-SCAPE" in e.value.args[0]
