# Configuration file for the Sphinx documentation builder.
# Full list of options: https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
import datetime

# -- Path setup ----------------------------------------------------------------

sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('../../'))

# -- Project information -------------------------------------------------------

project = 'MESH-Scripts-PyLib'
copyright = f'{datetime.datetime.now().year}, Fuad Yassin'
author = 'Fuad Yassin'
release = "1.0"

# -- General configuration -----------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',                  # Auto-document from docstrings
    'sphinx.ext.autosummary',              # Generate summary tables
    'sphinx.ext.viewcode',                 # Add links to source code
    'sphinx.ext.napoleon',                 # Support for NumPy/Google style docstrings
    'sphinx.ext.mathjax',                  # LaTeX math rendering
    'sphinx.ext.inheritance_diagram',      # Diagrams for class inheritance
    'sphinx.ext.graphviz',                 # Graphviz support
    'sphinx_automodapi.automodapi',        # Cleaner autodoc output for modules
    'sphinx_automodapi.smart_resolver',
    'sphinx_rtd_theme',                    # Read the Docs theme
    'nbsphinx',                            # Embed Jupyter notebooks
    'nbsphinx_link'                        # Link notebooks outside doc folder
]

autosummary_generate = True  # Auto-generate .rst files for documented modules/classes/functions

autodoc_default_options = {
    'members': True,
    'undoc-members': False,
    'private-members': False,
    'special-members': False,
    'inherited-members': False,
    'show-inheritance': True
}

# -- Napoleon config (for docstring style) -------------------------------------

napoleon_numpy_docstring = True
napoleon_google_docstring = False  # Set True if you're using Google style instead

# -- Jupyter notebook settings -------------------------------------------------

nbsphinx_execute = 'never'  # Don't run notebooks on doc build

# -- Template and build settings -----------------------------------------------

templates_path = ['_templates']
exclude_patterns = ['_build', '**.ipynb_checkpoints']

# -- HTML output options -------------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
pygments_style = 'sphinx'

html_context = {
    'display_github': True,
    'github_user': 'fuadyassin',
    'github_repo': 'MESH-Scripts-PyLib',
    'github_version': 'main/docs/source/'
}

# -- Inheritance diagram appearance --------------------------------------------

inheritance_graph_attrs = dict(rankdir="LR", size='"6.0, 8.0"', fontsize=12)
inheritance_node_attrs = dict(shape="ellipse", fontsize=12, color="blue", style="filled", fillcolor="lightgray")
