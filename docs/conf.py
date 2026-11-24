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
# The look follows the FESTIM brand: festim.css restyles furo from one block of
# tokens, mechanisms.css styles the drawings and animations of the theory page.
html_static_path = ["_static"]
html_css_files = ["festim.css", "mechanisms.css"]
html_favicon = "_static/festim-symbol.svg"
html_theme_options = {
    "light_logo": "festim-symbol.svg",
    "dark_logo": "festim-symbol-dark.svg",
}
# the Lucide "link" icon, shown next to a heading on hover
html_permalinks_icon = (
    '<svg aria-hidden="true" xmlns="http://www.w3.org/2000/svg" width="20" height="20"'
    ' viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"'
    ' stroke-linecap="round" stroke-linejoin="round">'
    '<path d="M10 13a5 5 0 0 0 7.54.54l3-3a5 5 0 0 0-7.07-7.07l-1.72 1.71"/>'
    '<path d="M14 11a5 5 0 0 0-7.54-.54l-3 3a5 5 0 0 0 7.07 7.07l1.71-1.71"/>'
    "</svg>"
)
exclude_patterns = ["_build"]
