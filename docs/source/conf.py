# Configuration file for the Sphinx documentation builder.
# Full list of options: https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
import datetime

# -- Path setup ----------------------------------------------------------------
# -- Path setup ----------------------------------------------------------------
# Add current and parent directories to the system path
# This allows Sphinx to locate your source modules (e.g., MESHpyPostProcessing)
sys.path.insert(0, os.path.abspath('../../src'))
#sys.path.insert(0, os.path.abspath('.'))
#sys.path.insert(0, os.path.abspath('../../'))

# -- Project information -------------------------------------------------------
project = 'MESH-Scripts-PyLib'
copyright = f'{datetime.datetime.now().year}, Fuad Yassin'
author = 'Fuad Yassin'
release = "1.0"   # Version of the project

# -- General configuration -----------------------------------------------------
# List of Sphinx extensions to load
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


# Automatically generate autosummary `.rst` files
# Auto-generate .rst files for documented modules/classes/functions
autosummary_generate = True  

# Settings for autodoc behavior
autodoc_default_options = {
    'members': True,                 # Include class/function members
    'undoc-members': False,          # Skip members without docstrings
    'private-members': False,        # Skip private members (with underscore prefix)
    'special-members': False,        # Skip special methods like __init__, __str__, etc.
    'inherited-members': False,      # Don’t include inherited members by default
    'show-inheritance': True         # Show class inheritance in docs
}

# -- Napoleon configuration (for docstring formatting) -------------------------
napoleon_numpy_docstring = True    # Enable NumPy-style docstrings
napoleon_google_docstring = False  # Set to True if using Google-style docstrings instead

# -- Jupyter notebook execution control ----------------------------------------
nbsphinx_execute = 'never'        # Don't execute notebooks during doc build

# -- Templates and source control ----------------------------------------------
templates_path = ['_templates']                   # Folder for custom templates
exclude_patterns = ['_build', '**.ipynb_checkpoints']  # Ignore build folders and notebook cache

# -- HTML output configuration -------------------------------------------------
html_theme = 'sphinx_rtd_theme'    # Use the popular Read the Docs theme
html_static_path = ['_static']     # Folder for custom CSS/JS (optional)
pygments_style = 'sphinx'          # Code highlighting style

# GitHub context for linking docs to repository
html_context = {
    'display_github': True,                   # Show GitHub link in top-right corner
    'github_user': 'fuadyassin',              # GitHub username
    'github_repo': 'MESH-Scripts-PyLib',      # Repository name
    'github_version': 'main/docs/source/'     # Branch + docs folder path
}

# -- Appearance of inheritance diagrams ----------------------------------------
# These settings customize Graphviz diagrams used by `sphinx.ext.inheritance_diagram`
inheritance_graph_attrs = dict(
    rankdir="LR",              # Left-to-right layout
    size='"6.0, 8.0"',         # Size of the output diagram
    fontsize=12
)
inheritance_node_attrs = dict(
    shape="ellipse",           # Node shape
    fontsize=12,
    color="blue",
    style="filled",
    fillcolor="lightgray"
)