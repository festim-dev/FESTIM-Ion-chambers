API reference
=============

Everything below is available from the top-level package, ``import ion_chambers as ic``,
except where a module path is shown.

Problems
--------

.. autoclass:: ion_chambers.ElectrostaticProblem
   :members: potential, initialise

.. autoclass:: ion_chambers.FixedPotentialBC

.. autoclass:: ion_chambers.PotentialExport

.. autoclass:: ion_chambers.ChargeTransportProblem
   :members: initialise

Geometry
--------

.. autoclass:: ion_chambers.Geometry
   :members:

.. autoclass:: ion_chambers.CoaxialChamber
   :members: from_volume, from_aspect_ratio, volume, create_mesh,
      is_anode, is_cathode

Gas
---

.. autoclass:: ion_chambers.Gas
   :members: number_density, thermal_voltage

.. autoclass:: ion_chambers.ChargeCarrier

.. autoclass:: ion_chambers.CarrierReaction
   :members: is_recombination

.. autodata:: ion_chambers.gas.Coefficient

Sources
-------

.. autoclass:: ion_chambers.IonisationSource
   :members:

.. autoclass:: ion_chambers.TritiumInGas

.. autofunction:: ion_chambers.source.lumped_projection

Exports
-------

.. autoclass:: ion_chambers.ChamberQuantity
   :members: compute

.. autoclass:: ion_chambers.ChamberSurfaceQuantity

.. autoclass:: ion_chambers.ChamberVolumeQuantity

.. autoclass:: ion_chambers.ChamberCurrent

.. autoclass:: ion_chambers.ElectrodeCurrent

.. autoclass:: ion_chambers.IdealCurrent

.. autoclass:: ion_chambers.WallEfficiency

.. autoclass:: ion_chambers.CollectionEfficiency

.. autoclass:: ion_chambers.PairGeneration

.. autoclass:: ion_chambers.Recombination

.. autoclass:: ion_chambers.SpaceChargeField

Stabilisation and FESTIM extensions
-----------------------------------

.. autoclass:: ion_chambers.DiffusionFloor
   :members:

.. autofunction:: ion_chambers.charge_transport_problem.cell_field_length

.. autoclass:: ion_chambers.AbsorbingElectrodeBC
   :members:

.. autoclass:: ion_chambers.MobilityElectromigration

.. autoclass:: ion_chambers.CellwiseDiffusionMaterial

Beta spectrum and wall losses
-----------------------------

.. automodule:: ion_chambers.beta
   :members:


Constants
---------

.. automodule:: ion_chambers.constants
   :members:
