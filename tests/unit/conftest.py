import os
import shutil
import tempfile
from . import DATA_DIR


nplinker_root_dir = os.path.join(tempfile.gettempdir(), "nplinker_unit_test")


def pytest_sessionstart(session):
    """Pytest hook to run before the entire test session starts.

    This hook makes sure the temporary directory `nplinker_root_dir` is created before any test
    starts. When running tests in parallel, the creation operation is done by the master process,
    and worker processes are not allowed to do it.

    For more about this hook, see:
    1. https://docs.pytest.org/en/stable/reference.html#_pytest.hookspec.pytest_sessionstart
    2. https://github.com/pytest-dev/pytest-xdist/issues/271#issuecomment-826396320
    """
    workerinput = getattr(session.config, "workerinput", None)
    # It's master process or not running in parallell when `workerinput` is None.
    if workerinput is None:
        if os.path.exists(nplinker_root_dir):
            shutil.rmtree(nplinker_root_dir)
        os.mkdir(nplinker_root_dir)
    # NPLinker setting `root_dir` must be a path that exists, so setting it to a temporary directory.
    os.environ["NPLINKER_ROOT_DIR"] = nplinker_root_dir
    # # Specify the config file via environment variable before importing nplinker in any test.
    os.environ["NPLINKER_CONFIG_FILE"] = str(DATA_DIR / "nplinker_local_mode.toml")


def pytest_sessionfinish(session):
    """Pytest hook to run after the entire test session finishes.

    This hook makes sure that temporary directory `nplinker_root_dir` is only removed after all
    tests finish. When running tests in parallel, the deletion operation is done by the master
    processs, and worker processes are not allowed to do it.

    """
    workerinput = getattr(session.config, "workerinput", None)
    if workerinput is None:
        shutil.rmtree(nplinker_root_dir)
