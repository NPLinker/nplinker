
???- Note "Requirements"
    - Linux, MacOS, or [WSL on Windows](https://learn.microsoft.com/en-us/windows/wsl/about)
        - For Windows without WSL enabled, please use NPLinker [docker image](https://hub.docker.com/r/nlesc/nplinker)
    - Python version ≥3.9


NPLinker is a python package that has both pypi packages and non-pypi packages as dependencies.
Install `nplinker` package as following:


```bash title="Install nplinker package"
# Check python version (≥3.9)
python --version

# Create a new virtual environment
python -m venv env          # (1)!
source env/bin/activate

# install nplinker package
pip install nplinker

# install nplinker non-pypi dependencies and databases
install-nplinker-deps
```

1. A virtual environment is ***required*** to install the the non-pypi dependencies. You can also use `conda` to create a new environment. But NPLinker is not available on conda yet.

## Install from source code

You can also install NPLinker from source code:

```bash title="Install from latest source code"
pip install git+https://github.com/nplinker/nplinker@dev  # (1)!
install-nplinker-deps
```

1. The `@dev` is the branch name. You can replace it with the branch name, commit or tag.
