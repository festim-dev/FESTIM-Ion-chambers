# Roadmap

The aim is a model for exploring chamber design: the linear range, memory effects and
the detection threshold. It has not been validated against measurements.

## Done

- Prototype model of a coaxial chamber in air, with response, shape and bias studies
  ([reference case](reference_case.md)).
- Package restructured along the lines of FESTIM: an electrostatic problem and a
  charge transport problem, with FESTIM-style subdomains, boundary conditions and
  exports, a test suite and this documentation.
- Verification of the transport problem: manufactured solutions on cartesian and
  cylindrical meshes, Boag's parallel-plate formula, and a transient run
  ([checks](reference_case.md#checks)).

## Next

1. **Carrier gases.** Worked examples for N₂ and Ar with electron transport that depends
   on E/N, and impurities that capture electrons. The model already accepts such
   coefficients; what is missing is vetted data and a choice of source for it.
2. **Wall losses from a mesh.** A mesh and its tags can already be read from a file
   and used in both problems. What is missing is a wall-loss Monte Carlo that traces
   rays against the mesh boundary instead of analytic cylinders.
3. **Space-charge coupling.** A class that iterates between the transport problem and
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
  mesh comes out too low (about 6% at 30 radial and 18 axial cells for the reference chamber).
  {class}`~ion_chambers.AbsorbingElectrodeBC` decides what crosses an electrode from
  the force on the carrier and does not have this error. It is opt-in for now. Before
  it becomes the default: a quantity that reports the charge leaving through the
  repelling electrode, so that this kind of loss is visible in any run.

## Tritium in the walls

Uptake, diffusion and release of tritium in the chamber wall, modelled with FESTIM as
a second problem. This gives the memory effect and the background drift after an
exposure. The two problems do not need to be solved together: charges cross the
chamber in milliseconds, while hydrogen moves through the wall over hours or more. The
wall problem is stepped in time and the chamber is solved as a steady state whenever a
reading is wanted. The [walkthrough](walkthrough.md#memory-effects-coupling-festim-to-the-chamber)
runs this chain for the tritium released into the gas. A full coupling adds:

1. **A source for the direct betas.** A `TritiumOnWall` source, a child of
   `IonisationSource`, takes the near-surface inventory from the FESTIM problem
   (Bq/m² within the beta range of the surface) and puts the ion pairs in the few
   millimetres of gas next to that wall. It uses the reverse of the wall-loss Monte
   Carlo: betas start in the metal and the energy they bring into the gas is counted.
   A transport problem already sums several sources, so it sits next to
   `TritiumInGas`.
2. **Surface physics in the wall.** At room temperature the memory of a chamber is
   often set by the surface (oxide layer, adsorbed water, isotope exchange) more than
   by diffusion in the bulk. FESTIM's surface kinetics and trapping describe this.
3. **A wall problem for each material.** The can and the anode are two FESTIM
   problems, or two subdomains, each with its own area and source.
4. **A balance for the gas.** A closed chamber, or a slow purge, follows
   $V \, dA/dt = \lambda \, \varphi_\mathrm{out} \, S_\mathrm{wall} - Q \, A$,
   stepped together with the wall problem.
5. **A class that runs the sequence**: step the wall, and solve the chamber when a
   reading is wanted. Because the coupling is one way, this is a loop and not an
   iteration.

## Later

- **Chambers that need a 3D mesh**, such as a wire cage.
- **Gas flow** through inlet and outlet ports, and the response to a step in
  concentration.
- **Columnar recombination** along beta tracks, and scattering in the wall-loss Monte
  Carlo.

## Improvements to the code

None of these gives a wrong result today: where the code cannot do something, it
raises an error. They are what limits the cases the package can describe.

- **Chamber current with any number of electrodes.**
  {class}`~ion_chambers.ChamberCurrent` takes its weighting potential from the
  potential divided by the bias, which holds for two electrodes without space charge,
  and raises otherwise. A small Laplace solve, with 1 on the readout electrode and 0
  on the others, removes the restriction. It is needed for guard rings and for
  space-charge coupling.
- **Coefficients that depend on more than E/N.** A mobility or a rate is a constant or
  a function of the reduced field. Temperature, pressure and position are not
  available to it. This comes before the carrier gases above.
- **A gas that can be reused.** A gas holds FESTIM species, which belong to one
  problem, so a study builds a new gas for every mesh.
- **Other emitters and sources on the walls.** The beta spectrum and the wall-loss
  Monte Carlo are written for tritium in the gas. Another emitter, or tritium held on
  a wall, needs an emitter object and a wall source. This comes before the
  wall coupling above.
- **3D meshes and MPI.** The cell length along the field and the wall-loss projection
  are written with numpy for triangles on one rank. This comes before the wire cage.
- **Electrostatics.** Only fixed potentials and one permittivity: no dielectric
  parts, no floating electrode.
- **Names of the electrodes.** `chamber.anode` and `chamber.cathode` name the rod and
  the can, which is wrong when the bias is reversed.
- **Upstream to FESTIM.** A diffusion field per species, a mobility on the
  electromigration term, and derived quantities on cylindrical meshes would each
  remove a workaround here.

## Known limitations of the code

- The mesh must be 2D and made of triangles, and the model runs on one MPI rank.
- `ChargeTransportProblem` is a child of FESTIM's `HydrogenTransportProblem`. FESTIM
  plans to make its discontinuous formulation the only one, and this package will
  have to follow.
- With more than two electrodes, set an
  {class}`~ion_chambers.AbsorbingElectrodeBC` for each carrier on each electrode and
  read the currents with {class}`~ion_chambers.ElectrodeCurrent`.
- A current computed on a single mesh with the default electrode conditions is too
  low. Extrapolate from two meshes, as the examples do.
- The model has not been compared with measurements yet.
