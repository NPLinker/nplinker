
!!!- Note "Requirements"
    - Linux, MacOS or [Windows with WSL](https://learn.microsoft.com/en-us/windows/wsl/)
    - Python version ≥3.9


NPLinker is a python package that has both pypi packages and non-pypi packages as dependencies. It 
requires <span style="color:red;">**~4.5GB**</span> of disk space to install all the dependencies. 

Install `nplinker` package as following:

```bash title="Install nplinker package"
# Create a new virtual environment and activate it
conda create -n npl-3.10 python=3.10  # (1)!
conda activate npl-3.10 

# install nplinker package (requiring ~300MB of disk space)
pip install --pre nplinker # (2)!

# install nplinker non-pypi dependencies and databases (~4GB)
install-nplinker-deps # (3)!
```

1. A virtual environment is ***required*** to install the non-pypi dependencies. It's recommended to use `conda` to manage python virtual environments. Check [conda documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) for installation instructions.
2. NPLinker v2 is still under development and released as [pre-release](https://pypi.org/project/nplinker/#history). To install the pre-release, you need the `--pre` option. 
3. Use `which install-nplinker-deps` command to check if the command `install-nplinker-deps` is in the activated python virtual environment. If not, make sure you have activated the virtual environment and/or reinstall the `nplinker` package.

## Install from source code

You can also install NPLinker from source code:

```bash title="Install from latest source code"
pip install git+https://github.com/nplinker/nplinker@dev  # (1)!
install-nplinker-deps
```

1. The `@dev` is the branch name. You can replace it with the branch name, commit or tag.
