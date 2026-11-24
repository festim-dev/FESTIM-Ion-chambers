# Roadmap

The aim is to reproduce the existing measurements first, then use the validated model
to explore chamber design: the linear range, memory effects and the detection
threshold.

## Done

- Prototype model of a coaxial chamber in air, with response, shape and bias studies
  ([results](results.md)).
- Measurements digitised from the published figures
  ([experimental data](experimental_data.md)).
- Package restructured along the lines of FESTIM: an electrostatic problem and a
  charge transport problem, with FESTIM-style subdomains, boundary conditions and
  exports, a test suite and this documentation.
- Verification of the transport problem: manufactured solutions on cartesian and
  cylindrical meshes, Boag's parallel-plate formula, and a transient run
  ([checks](results.md#checks)).

## Next

1. **Carrier gases.** Worked examples for N₂ and Ar with electron transport that depends
   on E/N, and impurities that capture electrons. The model already accepts such
   coefficients; what is missing is vetted data and a choice of source for it.
2. **Wall losses from a mesh.** A mesh and its tags can already be read from a file
   and used in both problems. What is missing is a wall-loss Monte Carlo that traces
   rays against the mesh boundary instead of analytic cylinders.
3. **The real 1 L monitor.** Its geometry, then the first comparison with the N₂
   calibration data.
4. **Space-charge coupling.** A class that iterates between the transport problem and
   the electrostatic problem, as FESTIM couples heat transfer and hydrogen transport,
   to extend the model beyond a few tens of Ci/m³.

## Numerics

- **Stabilisation without added diffusion.** The [diffusion floor](guide/transport.md#stabilisation)
  keeps the drift stable by raising the diffusion coefficient to tens or hundreds of
  times its physical value. This limits the model to first order in the cell size and
  exaggerates every effect of diffusion near a wall. Upwinding in a discontinuous
  formulation stabilises the drift without diffusing across the field. FESTIM is
  moving to that formulation, so the change is best made there and inherited here.
- **The default electrode condition.** With the floor, the zero density that the
  default sets on a repelling electrode acts as a sink, and the current on a single
  mesh comes out too low (9.5% at 30 radial cells for the reference chamber).
  {class}`~ion_chambers.AbsorbingElectrodeBC` decides what crosses an electrode from
  the force on the carrier and does not have this error. It is opt-in for now. Before
  it becomes the default: a quantity that reports the charge leaving through the
  repelling electrode, so that this kind of loss is visible in any run.

## Later

- **Tritium in the walls.** Uptake, diffusion and release of tritium in the chamber
  wall, modelled with FESTIM as a second problem. The two problems do not need to be
  solved together: charges cross the chamber in milliseconds, while hydrogen moves
  through the wall over hours or more. The plan is to step the wall problem in time and
  solve the chamber as a steady state whenever a reading is wanted, with the wall
  inventory entering as an extra ionisation source. This gives the memory effect and
  the background drift after an exposure.
- **20 cc wire cage**, which needs a 3D mesh.
- **Gas flow** through inlet and outlet ports, and the response to a step in
  concentration.
- **Columnar recombination** along beta tracks, and scattering in the wall-loss Monte
  Carlo.

## Known limitations of the code

- The numerical helpers assume a 2D triangular mesh.
- The wall-loss projection and the stabilisation run on one MPI rank only.
- `ChargeTransportProblem` is a child of FESTIM's `HydrogenTransportProblem`. FESTIM
  plans to make its discontinuous formulation the only one, and this package will
  have to follow.
- The electrodes are read as two surfaces, one above and one below the mean
  potential. More than two electrodes are not supported.
