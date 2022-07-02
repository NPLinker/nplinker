# `nplinker` developer documentation

If you're looking for user documentation, go [here](README.md).

## Code editor
We use [Visual Studio Code (VS Code)](https://code.visualstudio.com/) as code editor.
The VS Code settings for this project can be found in [.vscode](.vscode).
The settings will be automatically loaded and applied when you open the project with VS Code.
See [the guide](https://code.visualstudio.com/docs/getstarted/settings) for more info about workspace settings of VS Code.


## Setup

```shell
# Create a virtual environment, e.g. with
python3 -m venv venv

# activate virtual environment
source venv/bin/activate

# make sure to have a recent version of pip and setuptools
python3 -m pip install --upgrade pip setuptools

# install development dependencies
pip install --no-cache-dir --editable .[dev]

# install non-pypi dependecies
install-nplinker-deps
```

Afterwards check that the install directory is present in the `PATH` environment variable.

## Running the tests

There are two ways to run tests.

The first way requires an activated virtual environment with the development tools installed:

```shell
pytest -v
```

The second is to use `tox`, which can be installed separately (e.g. with `pip install tox`), i.e. not necessarily inside the virtual environment you use for installing `nplinker`, but then builds the necessary virtual environments itself by simply running:

```shell
tox
```

Testing with `tox` allows for keeping the testing environment separate from your development environment.
The development environment will typically accumulate (old) packages during development that interfere with testing; this problem is avoided by testing with `tox`.

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


We use [prospector](https://pypi.org/project/prospector/) for linting, [isort](https://pycqa.github.io/isort/) to sort imports, [autoflake](https://github.com/PyCQA/autoflake) to remove unused imports, and [yapf](https://github.com/google/yapf) for formatting, i.e. fixing readability of your code style.

Running the linters and formatters requires an activated virtual environment with the development tools installed.

```shell
# linting
prospector --profile .prospector.yml

# check import style for the project
isort --check .

# check import style and show any proposed changes as a diff
isort --check --diff .

# sort imports for the project
isort .

# remove unused imports for the project
# WARNING: python keyword `pass` will also be removed automatically
autoflake --in-place --remove-all-unused-imports  -r .

# format python style for the project
yapf -r -i .

# format python styple for specific python file
yapf -i filename.py
```

**Note:** We have set linter and formatter in VS Code [settings](.vscode),
so if you're using VS Code, you can also use its shortcut to do linting, sorting and formatting.
Besides, docstring style is also set, you can use [autoDocString](https://marketplace.visualstudio.com/items?itemName=njpwerner.autodocstring) to automatically generate docstrings.

## Generating the API docs

```shell
cd docs
make html
```

The documentation will be in `docs/_build/html`

If you do not have `make` use

```shell
sphinx-build -b html docs docs/_build/html
```

To find undocumented Python objects run

```shell
cd docs
make coverage
cat _build/coverage/python.txt
```

To [test snippets](https://www.sphinx-doc.org/en/master/usage/extensions/doctest.html) in documentation run

```shell
cd docs
make doctest
```

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
