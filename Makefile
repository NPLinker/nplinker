# .PHONY is used to declare that the targets are not files
.PHONY: install-dev clean clean-build clean-pyc clean-test clean-doc release build update-version sync-webapp-readme build-docs deploy-docs

help:
	@echo "Available commands to 'make':"
	@echo "  install-dev   : do an editable install of the NPLinker package for development" 
	@echo "  clean         : remove all build, test, coverage and Python artifacts"
	@echo "  clean-build   : remove build artifacts"
	@echo "  clean-pyc     : remove Python cache file artifacts"
	@echo "  clean-test    : remove test and coverage artifacts"
	@echo "  clean-doc     : remove doc build artifacts"
	@echo "  build         : build package"
	@echo "  release       : upload package to pypi"
	@echo "  build-docs    : build documentation for local development"
	@echo "  deploy-docs   : deploy documentation to GitHub Pages"
	@echo "  update-version: update NPLinker version (e.g. make update-version CURRENT_VERSION=0.1.0 NEW_VERSION=0.2.0)"

install-dev:	
	pip install -e ".[dev]"
	install-nplinker-deps

clean: clean-build clean-pyc clean-test clean-doc

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

clean-doc:
	rm -rf site

build: clean
	python -m build
	ls -l dist

release: update
	python -m twine upload dist/*


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

sync-webapp-readme:
	mkdir -p docs/webapp
	curl -sSf https://raw.githubusercontent.com/NPLinker/nplinker-webapp/main/README.md -o docs/webapp/readme.md

build-docs: sync-webapp-readme
	mkdocs serve

deploy-docs: sync-webapp-readme
ifndef version
	$(error version is not set. Usage: make deploy-docs version=YOUR_VERSION)
endif
	mike deploy -p -u $(version) latest