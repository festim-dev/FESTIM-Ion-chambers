"""Sphinx configuration for the ion_chambers documentation."""

project = "ion_chambers"
author = "ion_chambers contributors"
copyright = "2026, ion_chambers contributors"

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
html_title = "ion_chambers"
html_extra_path = ["_extra"]  # the prototype report, copied as is
exclude_patterns = ["_build", "_extra"]
