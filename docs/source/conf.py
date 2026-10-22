# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'chloroscan'
copyright = '2024-2025, Yuhao Tong'
author = 'Yuhao Tong, Rob Turnbull'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx_rtd_theme",
    # "nbsphinx",
    # "myst_parser",
    "sphinx.ext.mathjax",
    "sphinx.ext.githubpages",
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.napoleon",
    "sphinx_copybutton",
    "sphinx.ext.graphviz",
    'sphinx.ext.autosectionlabel',
    # 'sphinxcontrib.bibtex',
    # 'snakedoc'
]

templates_path = ['_templates']
exclude_patterns = []

html_theme_options = {
    "navigation_depth": 2,
}

autosectionlabel_prefix_document = True  # optionally namespace labels by filename

# Intersphinx options
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "torch": ("https://pytorch.org/docs/stable/index.html", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    # "pandas": ("https://pandas.pydata.org/pandas-docs/stable/", None),
    # "scikit-learn": ("https://scikit-learn.org/stable/", None),
    # "matplotlib": ("https://matplotlib.org/stable/", None),
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
