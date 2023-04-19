# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'GWSignal'
copyright = '2023, C.K.,J.F.N.,T.A.'
author = 'C.K.,J.F.N.,T.A.'
release = '0.0.1'

import os
import sys
import pathlib
sys.path.insert(0, os.path.abspath('../..'))
sys.path.insert(0, pathlib.Path(__file__).parents[3].resolve().as_posix())
sys.path.insert(0, os.path.abspath('../../..'))
sys.path.insert(0, os.path.abspath('../../../gwsignal/core'))
sys.path.insert(0, os.path.abspath('../../../gwsignal/models'))
# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx_rtd_theme',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.mathjax',
    'numpydoc',
    'nbsphinx',
    'sphinx.ext.autosectionlabel',
    'sphinx_tabs.tabs',
    "sphinx.ext.viewcode",
    'sphinx.ext.doctest',
    'sphinx.ext.napoleon'
]

templates_path = ['_templates']
exclude_patterns = []
numpydoc_show_class_members = False
source_suffix = ['.rst', '.md', '.txt']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
pygments_style = 'sphinx'
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
htmlhelp_basename = 'gwsignaldoc'
