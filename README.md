## Badges

| [fair-software.eu](https://fair-software.eu/) recommendations | |
| :-- | :--  |
| (1/5) code repository              | [![github repo badge](https://img.shields.io/badge/github-nplinker-000.svg?color=blue)](https://github.com/NPLinker/nplinker) |
| (2/5) license                      | [![github license badge](https://img.shields.io/github/license/NPLinker/nplinker)](https://github.com/NPLinker/nplinker) |
| (3/5) community registry           | [![pypi badge](https://img.shields.io/pypi/v/nplinker.svg?color=blue)](https://pypi.python.org/project/nplinker/) |
| (4/5) citation                     | [![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.4680218.svg)](https://doi.org/10.5281/zenodo.4680218) |
| (5/5) checklist                    | [![OpenSSF Scorecard](https://api.scorecard.dev/projects/github.com/NPLinker/nplinker/badge)](https://scorecard.dev/viewer/?uri=github.com/NPLinker/nplinker) |
| overall                       | [![fair-software badge](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F-green)](https://fair-software.eu) |
| **Other best practices**           |  |
| Documentation                      | [![Static Badge](https://img.shields.io/badge/Docs_Available-lightgreen)](https://nplinker.github.io/nplinker) [🔗](https://nplinker.github.io/nplinker)|
| Build & Test                       | [![build](https://github.com/NPLinker/nplinker/actions/workflows/build.yml/badge.svg)](https://github.com/NPLinker/nplinker/actions/workflows/build.yml) |
| Static analysis                    | [![workflow scq badge](https://sonarcloud.io/api/project_badges/measure?project=NPLinker_nplinker&metric=alert_status)](https://sonarcloud.io/dashboard?id=NPLinker_nplinker) |
| Coverage                           | [![workflow scc badge](https://sonarcloud.io/api/project_badges/measure?project=NPLinker_nplinker&metric=coverage)](https://sonarcloud.io/dashboard?id=NPLinker_nplinker) |
| Citation data consistency          | [![cffconvert](https://github.com/NPLinker/nplinker/actions/workflows/cffconvert.yml/badge.svg)](https://github.com/NPLinker/nplinker/actions/workflows/cffconvert.yml) |


# Natural Products Linker (NPLinker)
NPLinker is a python framework for data mining microbial natural products by integrating genomics and metabolomics data.

Original paper: [Ranking microbial metabolomic and genomic links in the NPLinker framework using complementary scoring functions](https://doi.org/10.1371/journal.pcbi.1008920).

## Setup and usage

### Requirement
- Linux or MacOS
- Python version ≥3.9


### Installation
NPLinker is a python package, using both pypi packages and non-pypi packages as dependencies. It 
requires <span style="color:red;">**~4.5GB**</span> of disk space to install all the dependencies. 

```shell
# Check python version (requiring ≥3.9)
python --version

# Create a new virtual environment
python -m venv env
source env/bin/activate

# install from nplinker releases (requiring ~300MB of disk space)
pip install --pre nplinker

# or install the latest from source code
pip install git+https://github.com/nplinker/nplinker@dev 

# install nplinker non-pypi dependencies and databases (~4GB)
install-nplinker-deps
```
A virtual environment is *required* to install the the non-pypi dependencies. You can also use `conda`
to manage python environments.

### Testing

To run the tests, you need to clone this repo and install the development dependencies:

```shell
# Create a new virtual environment
python -m venv env
source env/bin/activate

# Clone the repository and install the development dependencies
git clone https://github.com/NPLinker/nplinker.git
cd nplinker
pip install -e ".[dev]"
install-nplinker-deps
```

#### Unit tests

To run the unit tests, you can use the following command:

```shell
pytest
```
Pytest will use all available CPU cores to run the unit tests in parallel.

#### Integration tests

To run the integration tests, you can use the following command:

```shell
pytest -n1 tests/integration
```
The `-n1` is to use one CPU core to run the tests. Change it to `-n2` if you want to use two CPU cores to run in parallel.

### Usage

See the [documentation](https://nplinker.github.io/nplinker) for more information about how to use NPLinker.

## Contributing

If you want to contribute to the development of nplinker, have a look at the [contribution guidelines](CONTRIBUTING.md) and [README for developers](README.dev.md).
