"""Sphinx configuration for the FESTIM - Ion chambers documentation."""

project = "FESTIM - Ion chambers"
author = "FESTIM - Ion chambers contributors"
copyright = "2026, FESTIM - Ion chambers contributors"

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
]

# FEniCSx and FESTIM only install through conda. Mocking them lets the documentation
# build from a plain pip environment (see .github/workflows/docs.yml).
autodoc_mock_imports = [
    "festim",
    "dolfinx",
    "ufl",
    "basix",
    "mpi4py",
    "petsc4py",
    "scifem",
]
autodoc_member_order = "bysource"
autodoc_typehints = "description"
napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_use_ivar = True

myst_enable_extensions = ["dollarmath", "amsmath", "colon_fence"]
myst_heading_anchors = 3

html_theme = "furo"
html_title = "FESTIM - Ion chambers"
html_static_path = ["_static"]  # the drawings and animations of the theory page
html_css_files = ["mechanisms.css"]
exclude_patterns = ["_build"]
