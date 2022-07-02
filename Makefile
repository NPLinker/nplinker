.PHONY: clean clean-build clean-pyc clean-test release build

help:
	@echo "clean - remove all build, test, coverage and Python artifacts"
	@echo "clean-build - remove build artifacts"
	@echo "clean-pyc - remove Python file artifacts"
	@echo "clean-test - remove test and coverage artifacts"
	@echo "release - upload package to pypi"
	@echo "build - build package"
	@echo "update - update pip, build, twine packages"

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
	find . -name '__pycache__' -exec rm -fr {} +
	find . -name '*_cache' -exec rm -fr {} +

clean-test:
	rm -f .coverage

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
