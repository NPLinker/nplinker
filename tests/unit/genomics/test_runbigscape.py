from nplinker.genomics.bigscape.runbigscape import run_bigscape
from .. import DATA_DIR


class TestRunBigscape:
    def test_run_bigscape_v1(self, tmp_path):
        result = run_bigscape(
            antismash_path=DATA_DIR,
            output_path=tmp_path,
            extra_params="--help",
            version=1,
        )

        assert result is True

    def test_run_bigscape_v2(self, tmp_path):
        result = run_bigscape(
            antismash_path=DATA_DIR,
            output_path=tmp_path,
            extra_params="--help",
            version=2,
        )

        assert result is True
