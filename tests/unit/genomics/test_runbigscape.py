import os
import subprocess
import pytest
from nplinker.genomics import bigscape
from nplinker.genomics.bigscape import runbigscape
from .. import DATA_DIR


@pytest.mark.parametrize("version", ["1", "2"])
def test_run_bigscape(tmp_path, version):
    """Test whether BiG-SCAPE runs at all using the --help command"""
    result = bigscape.run_bigscape(
        antismash_path=tmp_path,
        output_path=tmp_path,
        extra_params="--help",
        version=version,
    )

    assert result is True


@pytest.mark.skipif(
    os.getenv("GITHUB_ACTIONS") == "true", reason="The test is time-consuming on CI"
)
@pytest.mark.parametrize("version", ["1", "2"])
def test_run_bigscape_small_dataset(tmp_path, version):
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
            version="3",
        )

    assert "version" in e.value.args[0]


@pytest.mark.parametrize("version", ["1", "2"])
def test_input_path_not_exist(tmp_path, version):
    with pytest.raises(FileNotFoundError) as e:
        bigscape.run_bigscape(
            antismash_path=tmp_path / "not_exist",
            output_path=tmp_path,
            extra_params="",
            version=version,
        )

    assert "antismash_path" in e.value.args[0]


@pytest.mark.parametrize("version", ["1", "2"])
def test_bad_parameters(tmp_path, version):
    with pytest.raises(RuntimeError) as e:
        bigscape.run_bigscape(
            antismash_path=tmp_path,
            output_path=tmp_path,
            extra_params="--this-is-not-a-real-argument",
            version=version,
        )

    assert "BiG-SCAPE" in e.value.args[0]


def test_v2_converts_underscores_to_hyphens(tmp_path, monkeypatch):
    """Test that underscores in option names are converted to hyphens for BiG-SCAPE v2."""
    calls = []

    def fake_run(args, **kwargs):
        calls.append(args)
        return subprocess.CompletedProcess(args, returncode=0)

    monkeypatch.setattr(runbigscape.subprocess, "run", fake_run)

    bigscape.run_bigscape(
        antismash_path=tmp_path,
        output_path=tmp_path,
        extra_params="--mibig_version 3.1 --include_singletons --gcf_cutoffs 0.30",
        version="2",
    )

    # Second call is the actual BiG-SCAPE run (first is the -h check)
    actual_args = calls[1]
    assert "--mibig-version" in actual_args
    assert "--include-singletons" in actual_args
    assert "--gcf-cutoffs" in actual_args
    assert "--mibig_version" not in actual_args


def test_v1_preserves_underscores(tmp_path, monkeypatch):
    """Test that underscores in option names are preserved for BiG-SCAPE v1."""
    calls = []

    def fake_run(args, **kwargs):
        calls.append(args)
        return subprocess.CompletedProcess(args, returncode=0)

    monkeypatch.setattr(runbigscape.subprocess, "run", fake_run)

    bigscape.run_bigscape(
        antismash_path=tmp_path,
        output_path=tmp_path,
        extra_params="--mibig_version 3.1",
        version="1",
    )

    actual_args = calls[1]
    assert "--mibig_version" in actual_args
