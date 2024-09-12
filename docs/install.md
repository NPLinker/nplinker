
???- Note "Requirements"
    - Linux, MacOS and [Windows with WSL](https://learn.microsoft.com/en-us/windows/wsl/)
    - Python version ≥3.9


NPLinker is a python package that has both pypi packages and non-pypi packages as dependencies. It 
requires <span style="color:red;">**~4.5GB**</span> of disk space to install all the dependencies. 

Install `nplinker` package as following:


```bash title="Install nplinker package"
# Check python version (≥3.9)
python --version

# Create a new virtual environment
python -m venv env          # (1)!
source env/bin/activate     # (2)! 

# install nplinker package (requiring ~300MB of disk space)
pip install --pre nplinker # (3)!

# install nplinker non-pypi dependencies and databases (~4GB)
install-nplinker-deps
```

1. A virtual environment is ***required*** to install the the non-pypi dependencies. You can also use `conda` to create a new environment. But NPLinker is not available on conda yet.
2. Check `pip` command and make sure it is provided by the activated virtual environment. 
3. NPLinker v2 is still under development and released as [pre-release](https://pypi.org/project/nplinker/#history). To install the pre-release, you need the `--pre` option. 

## Install from source code

You can also install NPLinker from source code:

```bash title="Install from latest source code"
pip install git+https://github.com/nplinker/nplinker@dev  # (1)!
install-nplinker-deps
```

1. The `@dev` is the branch name. You can replace it with the branch name, commit or tag.
