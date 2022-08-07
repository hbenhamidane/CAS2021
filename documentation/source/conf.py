# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup
import os
import pathlib
import sys
import sphinx_rtd_theme

sys.path.insert(0, pathlib.Path(__file__).parents[2].resolve().as_posix())
# sys.path.insert(0, os.path.abspath(os.path.join('..', '..')))


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Project Water'
copyright = '2022, Hisham Ben Hamidane, Ludovic Le Reste'
author = 'Hisham Ben Hamidane, Ludovic Le Reste'
release = '0.2'


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.autosummary',
              'sphinx_rtd_theme',
              'sphinx.ext.napoleon']

templates_path = ['_templates']
exclude_patterns = []

autosummary_generate = True
napoleon_google_docstring = False



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
