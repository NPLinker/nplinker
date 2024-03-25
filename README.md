## Badges

| [fair-software.eu](https://fair-software.eu/) recommendations | |
| :-- | :--  |
| (1/5) code repository              | [![github repo badge](https://img.shields.io/badge/github-nplinker-000.svg?color=blue)](https://github.com/NPLinker/nplinker) |
| (2/5) license                      | [![github license badge](https://img.shields.io/github/license/NPLinker/nplinker)](https://github.com/NPLinker/nplinker) |
| (3/5) community registry           | [![pypi badge](https://img.shields.io/pypi/v/nplinker.svg?color=blue)](https://pypi.python.org/project/nplinker/) ![Docker Image Version (latest by date)](https://img.shields.io/docker/v/nlesc/nplinker?arch=amd64&label=docker&sort=date)|
| (4/5) citation                     | [![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.4680218.svg)](https://doi.org/10.5281/zenodo.4680218) |
| (5/5) checklist                    | ![Static Badge](https://img.shields.io/badge/Coming_Soon-lightgrey) |
| how FAIR is                          | [![fair-software badge](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B-yellow)](https://fair-software.eu) |
| **Other best practices**           |  |
| Documentation                      | [![Static Badge](https://img.shields.io/badge/Docs_Available-lightgreen)](https://nplinker.github.io/nplinker) [🔗](https://nplinker.github.io/nplinker)|
| Build & Test                       | [![build](https://github.com/NPLinker/nplinker/actions/workflows/build.yml/badge.svg)](https://github.com/NPLinker/nplinker/actions/workflows/build.yml) |
| Static analysis                    | [![workflow scq badge](https://sonarcloud.io/api/project_badges/measure?project=NPLinker_nplinker&metric=alert_status)](https://sonarcloud.io/dashboard?id=NPLinker_nplinker) |
| Coverage                           | [![workflow scc badge](https://sonarcloud.io/api/project_badges/measure?project=NPLinker_nplinker&metric=coverage)](https://sonarcloud.io/dashboard?id=NPLinker_nplinker) |
| Citation data consistency          | [![cffconvert](https://github.com/NPLinker/nplinker/actions/workflows/cffconvert.yml/badge.svg)](https://github.com/NPLinker/nplinker/actions/workflows/cffconvert.yml) |
| **Downloads** | ![Docker Pulls](https://img.shields.io/docker/pulls/nlesc/nplinker?color=green&label=docker%20pulls) |


# Natural Products Linker (NPLinker)
NPLinker is a python framework for data mining microbial natural products by integrating genomics and metabolomics data.

Original paper: [Ranking microbial metabolomic and genomic links in the NPLinker framework using complementary scoring functions](https://doi.org/10.1371/journal.pcbi.1008920).

## Setup and usage

### Requirement
- Linux, MacOS, or Windows with WSL ([Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/about)) enabled
  - For Windows without WSL, please use the [docker image](https://hub.docker.com/r/nlesc/nplinker)
- Python version ≥3.9


### Installation
NPLinker is a Python package, you can install it as following:
```shell
# create a new virtual environment
python -m venv env
source env/bin/activate

# install nplinker package
pip install nplinker

# install nplinker non-pypi dependencies and databases
install-nplinker-deps
```

A virtual environment is *required* to install the the non-pypi dependencies. You can also use `conda`
to manage python environments.

### Usage

See the [documentation](https://nplinker.github.io/nplinker) for more information about how to use NPLinker.

## Contributing

If you want to contribute to the development of nplinker, have a look at the [contribution guidelines](CONTRIBUTING.md) and [README for developers](README.dev.md).
