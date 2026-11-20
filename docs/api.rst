API reference
=============

Everything below is available from the top-level package: ``import ion_chambers as ic``.

Model
-----

.. autoclass:: ion_chambers.IonChamberModel
   :members: initialise, run, concentrations

.. autoclass:: ion_chambers.Reading
   :members: as_dict

Geometry
--------

.. autoclass:: ion_chambers.Geometry
   :members:

.. autoclass:: ion_chambers.CoaxialChamber
   :members: from_volume, from_aspect_ratio, volume, default_axial_cells, build_mesh

Gas
---

.. autoclass:: ion_chambers.Gas
   :members: carrier, number_density, thermal_voltage

.. autoclass:: ion_chambers.ChargeCarrier

.. autoclass:: ion_chambers.Reaction
   :members: is_recombination

.. autodata:: ion_chambers.gas.Coefficient

Sources
-------

.. autoclass:: ion_chambers.IonisationSource
   :members:

.. autoclass:: ion_chambers.TritiumInGas

Beta spectrum and wall losses
-----------------------------

.. automodule:: ion_chambers.beta
   :members:

.. automodule:: ion_chambers.wall_loss
   :members:

Electric field
--------------

.. autoclass:: ion_chambers.PotentialSolver
   :members:

.. autoclass:: ion_chambers.LaplaceField

.. autoclass:: ion_chambers.SpaceChargePerturbation
   :members:

Stabilisation
-------------

.. autoclass:: ion_chambers.DiffusionFloor
   :members:

.. autofunction:: ion_chambers.stabilisation.cell_field_length

Readout functions
-----------------

.. automodule:: ion_chambers.readout
   :members:
   :exclude-members: Reading

Parameter studies
-----------------

.. autofunction:: ion_chambers.richardson

.. autofunction:: ion_chambers.response_curve

.. autofunction:: ion_chambers.save_json

FESTIM extensions
-----------------

.. autoclass:: ion_chambers.MobilityElectromigration

.. autoclass:: ion_chambers.CellwiseDiffusionMaterial

Helpers
-------

.. automodule:: ion_chambers.integration
   :members:

.. automodule:: ion_chambers.projection
   :members:

Constants
---------

.. automodule:: ion_chambers.constants
   :members:
