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

Check that everything works (about 5 s):

```bash
python scripts/quick_check.py
```

It solves the reference chamber at three activities and prints the readings.
