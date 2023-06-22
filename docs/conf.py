# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
from datetime import datetime
import importlib.metadata
importlib.metadata.distribution("sphinx-argparse")
print(importlib.metadata.distribution("sphinx-argparse").version)

# Source code directory, relative to this file

sys.path.insert(0, os.path.abspath("../"))
# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Immunopepper'
copyright = f'{datetime.now():%Y}, BMI lab, ETHZ'
author = 'Laurie Prélot, Matthias Hüser, Jiayu Chen, Andre Kahles, Gunnar Rätsch'
release = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.doctest',
    'sphinx.ext.todo',
    'sphinx.ext.autosummary',
    'sphinxarg.ext',

]

templates_path = ['_templates']
exclude_patterns = []
source_suffix = ['.rst', '.md']


todo_include_todos = True

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

