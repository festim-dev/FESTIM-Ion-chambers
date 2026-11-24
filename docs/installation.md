# Installation

FESTIM and DOLFINx install through conda, not pip, so the package comes with a conda
environment file.

```bash
git clone https://github.com/festim-dev/ion_chambers
cd ion_chambers
conda env create -f environment.yml
conda activate ion-chambers
pip install -e .
```

This gives FESTIM 2.2 on DOLFINx 0.11 from conda-forge, plus the tools for testing
and building this documentation.

Check that everything works by running the tests (about 10 s):

```bash
python -m pytest
```
