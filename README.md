# Natural Products Linker (NPLinker)

NPLinker is a tool that aims to address the significant bottleneck that exists in the realization of the potential of genome-led metabolite discovery, namely the slow manual matching of predicted biosynthetic gene clusters (BGCs) with metabolites produced during bacterial culture; linking phenotype to genotype. NPLinker implements a new data-centric approach to alleviate this linking problem by searching for patterns of strain presence and absence between groups of similar spectra (molecular families; MF) and groups of similar BGCs (gene cluster families; GCF). Searches can be performed using a number of available analysis methods employed in isolation or together. 

## Getting Started - Web application

The simplest option for most users will be to run NPLinker as a Dockerised web application, removing the need to install Python and other dependencies as well as providing a graphical interface. The web application can also be used without Docker if required, but it's much simpler to rely on the Docker version. Assuming you have Docker installed already, see the [wiki page](https://github.com/sdrogers/nplinker/wiki/WebappInstallation) for detailed installation and usage instructions. 

## Getting Started - Python

NPLinker is implemented as a Python module. If you want maximum flexibility and are comfortable working in Python, you can set up a local environment as follows:
```
git clone https://github.com/sdrogers/nplinker
cd nplinker
pipenv sync
pipenv shell
```
See the example Jupyter notebook (notebooks/nplinker_demo1.ipynb) for a guided introduction to the NPLinker API which shows how to load and examine a dataset. 

API documentation is still being updated but the latest version is available at [readthedocs](https://nplinker.readthedocs.io/en/latest/).
