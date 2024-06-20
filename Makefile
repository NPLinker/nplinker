.PHONY: clean clean-build clean-pyc clean-test release build update-version

help:
	@echo "clean - remove all build, test, coverage and Python artifacts"
	@echo "clean-build - remove build artifacts"
	@echo "clean-pyc - remove Python file artifacts"
	@echo "clean-test - remove test and coverage artifacts"
	@echo "release - upload package to pypi"
	@echo "build - build package"
	@echo "update - update pip, build, twine packages"
	@echo "update-version - update NPLinker version. Usage: make update-version CURRENT_VERSION=0.1.0 NEW_VERSION=0.2.0"

clean: clean-build clean-pyc clean-test

clean-build:
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-pyc:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '*__pycache__' -exec rm -fr {} +
	find . -name '*_cache' -exec rm -fr {} +

clean-test:
	rm -f .coverage*
	rm -f coverage.xml

build: clean
	python -m build
	ls -l dist

release: update
	python -m twine upload dist/*

update:
	pip install --upgrade pip build twine

venv:
	python -m venv venv

clean-venv:
	rm -rf venv


# Define the files to update version
FILES := src/nplinker/__init__.py pyproject.toml CITATION.cff

# Rule to update the version in the specified files
update-version:
ifndef CURRENT_VERSION
	$(error CURRENT_VERSION is not provided. Usage: make update-version CURRENT_VERSION=0.1.0 NEW_VERSION=0.2.0)
endif
ifndef NEW_VERSION
	$(error NEW_VERSION is not provided. Usage: make update-version CURRENT_VERSION=0.1.0 NEW_VERSION=0.2.0)
endif
	@for file in $(FILES); do \
		if ! grep -qE "__version__ = \"$(CURRENT_VERSION)\"|version = \"$(CURRENT_VERSION)\"|version: \"$(CURRENT_VERSION)\"" $$file; then \
			echo "Error: Current version $(CURRENT_VERSION) not found in $$file"; \
			exit 1; \
		fi; \
	done

	@echo "Updating version from $(CURRENT_VERSION) to $(NEW_VERSION) for following files:"
	@for file in $(FILES); do \
		echo "  $$file"; \
		if [ "$(shell uname)" = "Darwin" ]; then \
			sed -i '' -e 's/__version__ = "$(CURRENT_VERSION)"/__version__ = "$(NEW_VERSION)"/' \
				-e 's/version = "$(CURRENT_VERSION)"/version = "$(NEW_VERSION)"/' \
				-e 's/version: "$(CURRENT_VERSION)"/version: "$(NEW_VERSION)"/' $$file; \
		else \
			sed -i'' -e 's/__version__ = "$(CURRENT_VERSION)"/__version__ = "$(NEW_VERSION)"/' \
				-e 's/version = "$(CURRENT_VERSION)"/version = "$(NEW_VERSION)"/' \
				-e 's/version: "$(CURRENT_VERSION)"/version: "$(NEW_VERSION)"/' $$file; \
		fi; \
	done
	@echo "Version update complete."