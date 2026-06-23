# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

sys.path.insert(0, os.path.abspath("../../ringtail/"))
sys.path.insert(0, os.path.abspath("../.."))  # repo root, so `import ringtail` resolves
from ringtail.ringtailoptions import ringtail_defaults


def _fmt_default(v):
    # empty list/tuple/string and None all render as "none" — an empty
    # ".. |x| replace::" is an invalid (warning-raising) substitution.
    if isinstance(v, (tuple, list)):
        return ", ".join(map(str, v)) if v else "none"
    s = "none" if v is None else str(v)
    return s or "none"


# turn every default into a substitution, e.g. |default_storage_type| -> duckdb
rst_prolog = "\n".join(
    f".. |default_{key}| replace:: {_fmt_default(val)}"
    for key, val in ringtail_defaults().items()
)


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "ringtail"
copyright = "2024, Forli lab"
author = "Forli lab"
_version_file = os.path.join(
    os.path.dirname(__file__), "..", "..", "ringtail", "_version.py"
)
_ns = {}
with open(_version_file) as _fh:
    exec(_fh.read(), _ns)
release = _ns["__version__"]
version = ".".join(release.split(".")[:2])

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx_design",
]

# Prefix autosection labels with the document name so identical section titles in
# different files (e.g. "Scoring filters" in api.rst and cmdline.rst) don't collide.
# All cross-references in the docs use explicit ".. _label:" targets, so this is safe.
autosectionlabel_prefix_document = True

# The changelog (changes.rst) intentionally repeats "Bug fixes"/"Enhancements" headings
# per version block; nothing references those auto-labels, so silence the duplicates.
suppress_warnings = ["autosectionlabel.changes"]

# Render docstring "Attributes:" sections as per-class :ivar: fields instead of global
# attribute targets, so attributes that share a name across classes (e.g. `cutoff` in the
# two clustermanager classes) don't create an ambiguous cross-reference.
napoleon_use_ivar = True

templates_path = ["_templates"]
exclude_patterns = []

pygments_style = "sphinx"


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "shibuya"
html_static_path = ["_static"]
html_css_files = ["custom.css"]

# Shibuya theme options (AutoDock-blue accent; light/dark toggle is built in).
html_theme_options = {
    "accent_color": "blue",
    "github_url": "https://github.com/forlilab/Ringtail",
}

autodoc_mock_imports = [
    "matplotlib",
    "meeko",
    "pandas",
    "rdkit",
    "numpy",
    "multiprocess",
]
