
!!!- Note "Requirements"
    - Linux, MacOS or [Windows with WSL](https://learn.microsoft.com/en-us/windows/wsl/)
    - Python version ≥3.11


NPLinker is a python package that has both pypi packages and non-pypi packages as dependencies. It 
requires <span style="color:red;">**~4.5GB**</span> of disk space to install all the dependencies. 

Install `nplinker` package as following:

```bash title="Install nplinker package"
# Create a new virtual environment and activate it
conda create -n npl-3.11 python=3.11  # (1)!
conda activate npl-3.11 

# install nplinker package (requiring ~300MB of disk space)
pip install nplinker

# install nplinker non-pypi dependencies and databases (~4GB)
install-nplinker-deps # (2)!
```

1. A virtual environment is ***required*** to install the non-pypi dependencies. It's recommended to use `conda` to manage python virtual environments. Check [conda documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) for installation instructions.
2. Use `which install-nplinker-deps` command to check if the command `install-nplinker-deps` is in the activated python virtual environment. If not, make sure you have activated the virtual environment and/or reinstall the `nplinker` package.

## Install from source code

You can also install NPLinker from source code:

```bash title="Install from latest source code"
pip install git+https://github.com/nplinker/nplinker@dev  # (1)!
install-nplinker-deps
```

1. The `@dev` is the branch name. You can replace it with the branch name, commit or tag.
