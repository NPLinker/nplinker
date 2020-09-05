# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('../prototype'))


# -- Project information -----------------------------------------------------

project = 'nplinker'
copyright = '2020, various'
author = 'various'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.coverage',
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.napoleon',
]

# (napoleon is an extension that allows use of non-RST docstrings)
# support Google style docstrings only
napoleon_google_docstring = True
napoleon_numpy_docstring = False
# list __init__ docstrings separately from class documentation
napoleon_include_init_with_doc = True
# don't include private member docstrings
napoleon_include_private_with_doc = False
# don't include __name__ members with docstrings in the documentation
napoleon_include_special_with_doc = False
# use :param: role for each func parameter if True, else single :parameters: role
napoleon_use_param = True
# use :keyword: role for each func keyword parameter if True, else single :keyword arguments role:
napoleon_use_keyword = True
# if True, format return type separately from description (vs inline)
napoleon_use_rtype = False

# include both class-level and __init__ method docstrings 
#autoclass_content = 'both'
# fake imports for autodoc
autodoc_mock_imports = ['Bio', 'toml', 'xdg', 'httpx', 'progress', 'sortedcontainers', 'seaborn', 'matplotlib', 'pandas']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

master_doc = 'index'


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
