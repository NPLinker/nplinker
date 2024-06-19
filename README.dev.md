# `nplinker` developer documentation

If you're looking for user documentation, go [here](README.md).

## Code editor
We use [Visual Studio Code (VS Code)](https://code.visualstudio.com/) as code editor.

The VS Code Profile for this project is [vscode/nplinker.code-profile](vscode/nplinker.code-profile), 
which contains the settings, extensions and snippets for the project. To use the profile, you must
first import it by clicking the following menus: `Code` -> `Settings` -> `Profiles` -> `Import Profile...`. 
Then select the file [vscode/nplinker.code-profile](vscode/nplinker.code-profile) to import the profile.
VS Code will take a while to install the extensions and apply the settings. Want more info? See 
[vscode profiles guide](https://code.visualstudio.com/docs/editor/profiles).


If you want to add more settings, you can update the workspace settings, see [the guide](https://code.visualstudio.com/docs/getstarted/settings) for more info.


## Setup

We use Python 3.10 for development environment.

```shell
# Create a virtual environment, e.g. with
python3 -m venv venv

# activate virtual environment
source venv/bin/activate

# make sure to have a recent version of pip and setuptools
python3 -m pip install --upgrade pip setuptools

# install development dependencies
pip install --no-cache-dir --editable ".[dev]"

# install non-pypi dependencies
install-nplinker-deps
```

Afterwards check that the install directory is present in the `PATH` environment variable.

You can also use [conda](https://docs.conda.io/projects/conda/en/stable/) to manage python environments.

## Running the tests

**Run unit tests with**
```shell
pytest
# or
pytest -n auto tests/unit
```
Parallel testing is supported with `pytest-xdist` plugin. To run tests in parallel, use the `-n`
option, e.g. `-n auto` to run tests in parallel with the number of CPUs available.

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
mypy path-to-source-code
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
The docs are published on github pages. We use [mike](https://github.com/jimporter/mike)
to deploy the docs to the `gh-pages` branch and to manage the versions of docs.

For example, to deploy the version 2.0 of the docs to the `gh-pages` branch and make it the latest
version, run:
```shell
mike deploy -p -u 2.0 latest
```
If you are not happy with the changes you can run `mike delete [version]`.
All these mike operations will be recorded as git commits of branch `gh-pages`.

 `mike serve` is used to check all versions committed to branch `gh-pages`, which is for checking
 the production website. If you have changes but not commit them yet, you should use `mkdocs serve`
 instead of  `mike serve` to check them.



## Versioning

Bumping the version across all files is done with [bumpversion](https://github.com/c4urself/bump2version), e.g.

```shell
bumpversion major
bumpversion minor
bumpversion patch

bumpversion --current-version 2.0.0 --new-version 2.1.0 fakepart
```

## Making a release

This section describes how to make a release in 3 parts:

1. preparation
1. making a release on PyPI
1. making a release on GitHub

### (1/3) Preparation

1. Update the <CHANGELOG.md> (don't forget to update links at bottom of page)
2. Verify that the information in `CITATION.cff` is correct, and that `.zenodo.json` contains equivalent data
3. Make sure the [version has been updated](#versioning).
4. Run the unit tests with `pytest -v`

### (2/3) PyPI

In a new terminal, without an activated virtual environment or an env directory:

```shell
# prepare a new directory
cd $(mktemp -d nplinker.XXXXXX)

# fresh git clone ensures the release has the state of origin/main branch
git clone https://github.com/NPLinker/nplinker .

# prepare a clean virtual environment and activate it
python3 -m venv env
source env/bin/activate

# make sure to have a recent version of pip and setuptools
python3 -m pip install --upgrade pip setuptools

# install runtime dependencies and publishing dependencies
python3 -m pip install --no-cache-dir .
python3 -m pip install --no-cache-dir .[publishing]

# clean up any previously generated artefacts
rm -rf nplinker.egg-info
rm -rf dist

# create the source distribution and the wheel
python3 -m build

# upload to test pypi instance (requires credentials)
twine upload --repository-url https://test.pypi.org/legacy/ dist/*
```

Visit
[https://test.pypi.org/project/nplinker](https://test.pypi.org/project/nplinker)
and verify that your package was uploaded successfully. Keep the terminal open, we'll need it later.

In a new terminal, without an activated virtual environment or an env directory:

```shell
cd $(mktemp -d nplinker-test.XXXXXX)

# prepare a clean virtual environment and activate it
python3 -m venv env
source env/bin/activate

# make sure to have a recent version of pip and setuptools
pip install --upgrade pip setuptools

# install from test pypi instance:
python3 -m pip -v install --no-cache-dir \
--index-url https://test.pypi.org/simple/ \
--extra-index-url https://pypi.org/simple nplinker
```

Check that the package works as it should when installed from pypitest.

Then upload to pypi.org with:

```shell
# Back to the first terminal,
# FINAL STEP: upload to PyPI (requires credentials)
twine upload dist/*
```

### (3/3) GitHub

Don't forget to also make a [release on GitHub](https://github.com/NPLinker/nplinker/releases/new). If your repository uses the GitHub-Zenodo integration this will also trigger Zenodo into making a snapshot of your repository and sticking a DOI on it.
