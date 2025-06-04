# `nplinker` developer documentation

If you're looking for user documentation, go [here](README.md).

## Code editor
We use [Visual Studio Code (VS Code)](https://code.visualstudio.com/) as code editor.

The VS Code Profile for this project is [vscode/nplinker.code-profile](vscode/nplinker.code-profile), which contains the settings, extensions and snippets for the project. 

To use the profile, you must first import it by clicking the following menus: `Code` -> `Settings` -> `Profiles` -> `Import Profile...`. 
Then select the file [vscode/nplinker.code-profile](vscode/nplinker.code-profile) to import the profile.
VS Code will take a while to install the extensions and apply the settings. Want more info? See [vscode profiles guide](https://code.visualstudio.com/docs/editor/profiles).

If you want to add more settings, you can update the workspace settings, see [the guide](https://code.visualstudio.com/docs/getstarted/settings) for more info.


## Setup

We use Python 3.10 for development environment.

```shell
# Create a virtual environment
conda create -n npl-dev python=3.10

# activate virtual environment
conda activate npl-dev

# Clone the repository
git clone https://github.com/NPLinker/nplinker.git
cd nplinker

# install development dependencies
pip install -e ".[dev]"

# install non-pypi dependencies
install-nplinker-deps
```


## Running the tests

**Run unit tests with**
```shell
pytest
# or
pytest -n 2 tests/unit
```
Parallel testing is supported with `pytest-xdist` plugin. To run tests in parallel, use the `-n` option, e.g. `-n 2` to run tests in parallel with 2 CPUs. 
By default, `pytest` will use all available CPUs to run the tests in parallel.

**Run integration tests with**
```shell
pytest -n 0 tests/integration
```
`-n 0` means no parallel testing.


### Test coverage

In addition to just running the tests to see if they pass, they can be used for coverage statistics, i.e. to determine how much of the package's code is actually executed during tests.
In an activated virtual environment with the development tools installed, inside the package directory, run:

```shell
coverage run
```

This runs tests and stores the result in a `.coverage` file.
To see the results on the command line, run

```shell
coverage report
```

`coverage` can also generate output in HTML and other formats; see `coverage help` for more information.


## Linting and formatting

We use [ruff](https://docs.astral.sh/ruff/) for linting, sorting imports and formatting code. The configurations of `ruff` are set in [pyproject.toml](pyproject.toml) file.

Running the linters and formatters requires an activated virtual environment with the development tools installed.

```shell
# Lint all files in the current directory.
ruff check .

# Lint all files in the current directory, and fix any fixable errors.
ruff check . --fix

# Format all files in the current directory
ruff format .

# Format a single python file
ruff format filename.py
```

## Static typing

We use [inline type annotation](https://typing.readthedocs.io/en/latest/source/libraries.html#how-to-provide-type-annotations) for static typing rather than stub files (i.e. `.pyi` files).

Since Python 3.10 is used as dev environment and NPLinker must support Python version ≥3.9, you may see various typing issues at runtime. Here is [a guide to solve the potential runtime issues](https://mypy.readthedocs.io/en/stable/runtime_troubles.html).

By default, we use `from __future__ import annotations` at module level to stop evaluating annotations at function definition time (see [PEP 563](https://peps.python.org/pep-0563/)), which would solve most of compatibility issues between different Python versions. Make sure you're aware of the [caveats](https://mypy.readthedocs.io/en/stable/runtime_troubles.html#future-annotations-import-pep-563).

We use [Mypy](http://mypy-lang.org/) as static type checker:

```
# install mypy
pip install mypy

# run mypy
mypy src/nplinker
```

Mypy configurations are set in [pyproject.toml](pyproject.toml) file.

For more info about static typing and mypy, see:
- [Static typing with Python](https://typing.readthedocs.io/en/latest/index.html#)
- [Mypy doc](https://mypy.readthedocs.io/en/stable/)


## Docs
We use [MkDocs](https://www.mkdocs.org/) and its theme [Material for MkDocs](https://squidfunk.github.io/mkdocs-material/)
to generate documentations. The configurations of MkDocs are set in [mkdocs.yml](mkdocs.yml) file.

To watch the changes of current doc in real time, run:
```shell
mkdocs serve
# or to watch src and docs directories
mkdocs serve -w docs -w src
```
Then open your browser and go to [http://127.0.0.1:8000/](http://127.0.0.1:8000/).

### Publishing the docs

Documentation is published to GitHub Pages using [mike](https://github.com/jimporter/mike), which also manages versioning on the `gh-pages` branch.

#### Deploying a new version

To deploy version `2.0` of the docs and mark it as the latest, run:

```bash
make deploy-docs version=2.0
```

This command does the following:

* Fetches the latest README from the [`nplinker-webapp`](https://github.com/NPLinker/nplinker-webapp) repository.
* Builds and deploys the documentation to the `gh-pages` branch.
* Updates the `latest` alias to point to this version.
* Creates a commit on the `gh-pages` branch to record the deployment.

If you want to undo a deployment, you can run:

```bash
mike delete 2.0
```

#### Previewing the docs

* To preview all committed versions from `gh-pages`, use:

  ```bash
  mike serve
  ```

* To preview your local, uncommitted changes, use:

  ```bash
  make build-docs
  ```

> `make build-docs` will also update the webapp README before serving the docs.

## Versioning

Updating the version of the NPLinker package is done with make command `update-version`, e.g.

```shell
make update-version CURRENT_VERSION=0.0.1 NEW_VERSION=0.0.2
```

This command will update the version in the following files:
- `src/nplinker/__init__.py`
- `pyproject.toml`
- `CITATION.cff`

## Making a release

This section describes how to make a release in 2 parts:

1. Create Github release
2. Publish to Pypi

### (1/2) Create Github release

We use the Github action [Draft or publish Github release
](https://github.com/NPLinker/nplinker/actions/workflows/publish_gh_release.yml) to create a Github release.

Click the right corner `Run workflow` button, then fill in the current version number and new version number, and choose `publish` to publish Github release, then click the `Run workflow` button. 

The action will first update the version with the command `make update-version`. Then it will generate a release notes and update the `CHANGELOG.md` file with the notes. After that, the action will commit and push the changes. In the end, the action will create a Github release with the new version number and create a tag for the release.

After the action is finished successfully, you can go to the [release page](https://github.com/NPLinker/nplinker/releases) to check the release.

This repository uses the GitHub-Zenodo integration, the new Github release will trigger Zenodo into making a snapshot of the repository and sticking a DOI on it. Check the [Zenodo page](https://zenodo.org/records/14723594) to see the new snapshot.

### (2/2) Publish to Pypi

You can publish the package to pypi with the following steps:

```shell
# Go to your local nplinker repository
cd path-to-nplinker-repo

# Clean the repository
make clean

# Build the source distribution and the wheel
make build

# Publish to pypi
make release
```

After publishing to pypi, you can check the [pypi page](https://pypi.org/project/nplinker/#history) to see the new version.