## Badges

| [fair-software.eu](https://fair-software.eu/) recommendations | |
| :-- | :--  |
| (1/5) code repository              | [![github repo badge](https://img.shields.io/badge/github-nplinker-000.svg?color=blue)](https://github.com/NPLinker/nplinker) |
| (2/5) license                      | [![github license badge](https://img.shields.io/github/license/NPLinker/nplinker)](https://github.com/NPLinker/nplinker) |
| (3/5) community registry           | [![pypi badge](https://img.shields.io/pypi/v/nplinker.svg?color=blue)](https://pypi.python.org/project/nplinker/) ![Docker Image Version (latest by date)](https://img.shields.io/docker/v/nlesc/nplinker?arch=amd64&label=docker&sort=date)|
| (4/5) citation                     | [![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.4680218.svg)](https://doi.org/10.5281/zenodo.4680218) |
| (5/5) checklist                    | Coming soon |
| how FAIR is                          | [![fair-software badge](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B-yellow)](https://fair-software.eu) |
| **Other best practices**           |  |
| Documentation                      | [![Documentation Status](https://readthedocs.org/projects/nplinker/badge/?version=latest)](https://nplinker.readthedocs.io/en/latest/?badge=latest) |
| Build & Test                       | [![build](https://github.com/NPLinker/nplinker/actions/workflows/build.yml/badge.svg)](https://github.com/NPLinker/nplinker/actions/workflows/build.yml) |
| Static analysis                    | [![workflow scq badge](https://sonarcloud.io/api/project_badges/measure?project=NPLinker_nplinker&metric=alert_status)](https://sonarcloud.io/dashboard?id=NPLinker_nplinker) |
| Coverage                           | [![workflow scc badge](https://sonarcloud.io/api/project_badges/measure?project=NPLinker_nplinker&metric=coverage)](https://sonarcloud.io/dashboard?id=NPLinker_nplinker) |
| Citation data consistency          | [![cffconvert](https://github.com/NPLinker/nplinker/actions/workflows/cffconvert.yml/badge.svg)](https://github.com/NPLinker/nplinker/actions/workflows/cffconvert.yml) |
| **Downloads** | ![Docker Pulls](https://img.shields.io/docker/pulls/nlesc/nplinker?color=green&label=docker%20pulls) |


# Natural Products Linker (NPLinker)

NPLinker aims to address the significant bottleneck that exists in the realization of the potential of genome-led metabolite discovery, namely the slow manual matching of predicted biosynthetic gene clusters (BGCs) with metabolites produced during bacterial culture; linking phenotype to genotype.

NPLinker implements a new data-centric approach to alleviate this linking problem by searching for patterns of strain presence and absence between groups of similar spectra (molecular families; MF) and groups of similar BGCs (gene cluster families; GCF). Searches can be performed using a number of available analysis methods employed in isolation or together.


Currently available analysis methods (scoring methods):
- **Metcalf** (standardised): see [Hjörleifsson Eldjárn G, et al. (2021)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008920), [Doroghazi JR, et al. (2014)](https://www.nature.com/articles/nchembio.1659), and the [demo](notebooks/nplinker_demo1.ipynb).
- **Rosetta**: see [Hjörleifsson Eldjárn G, et al. (2021)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008920) and [Soldatou S, et al. (2021)](https://www.mdpi.com/1660-3397/19/2/103).
- **NPClassScore**: see the preprint [Louwen JJR, et al. (2022)](https://www.researchsquare.com/article/rs-1391827/v2), and the [demo](notebooks/npclassscore_linking/NPClassScore_demo.ipynb).

## Setup and usage

NPLinker is a Python package, you can install it as following:
```
# create a new virtual environment
python -m venv env
source env/bin/activate

# install nplinker package
pip install nplinker

# install nplinker non-pypi dependencies and databases
install-nplinker-deps
```

Due to hardware requirements of some non-pypi dependecies:
- NPLinker can only be installed on `Linux` and `MacOS (Intel chip)`
- `MacOS(Apple Silicon, M1/M2 chip)` user should execute the install commands in a [Rosseta-enabled terminal](https://support.apple.com/en-us/HT211861).
- For `Windows` users, please use our [docker image](https://hub.docker.com/r/nlesc/nplinker).


See the example in [Jupyter notebook](notebooks/nplinker_demo1.ipynb) for a guided introduction to the NPLinker API which shows how to load and examine a dataset. Other notebooks are present showcasing other scoring methods, like for [NPClassScore](notebooks/npclassscore_linking/NPClassScore_demo.ipynb).

If you want to visualize and manipulate NPLinker predictions, check [NPLinker Webapp](https://github.com/NPLinker/nplinker-webapp) for more info.

## Contributing

If you want to contribute to the development of nplinker, have a look at the [contribution guidelines](CONTRIBUTING.md) and [README for developers](README.dev.md).
