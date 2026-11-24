# Installation

FESTIM and DOLFINx install through conda, not pip, so the package comes with a conda
environment file.

```bash
git clone https://github.com/festim-dev/FESTIM-Ion-chambers
cd FESTIM-Ion-chambers
conda env create -f environment.yml
conda activate ion-chambers-env
pip install -e .
```

This gives FESTIM 2.2 on DOLFINx 0.11 from conda-forge, plus the tools for testing
and building this documentation.

Check that everything works by running the tests:

```bash
python -m pytest
```
