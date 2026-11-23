"""Integration tests: full solves of small problems."""

from mpi4py import MPI

import dolfinx
import festim as F
import io4dolfinx
import numpy as np
import pytest
from conftest import ACTIVITIES
from prototype_case import air, reference_chamber, solve_potential

import ion_chambers as ic

pytestmark = pytest.mark.slow

MICROCURIE = ACTIVITIES.index(3.7e4)  # 1 uCi/m3
TEN_MICROCURIE = ACTIVITIES.index(3.7e5)
CURIE = ACTIVITIES.index(3.7e10)  # 1 Ci/m3
THOUSAND_CURIE = ACTIVITIES.index(3.7e13)


def electron_ion_gas(mobility=0.04) -> ic.Gas:
    """A non-attaching gas: electrons and positive ions only."""
    electrons = ic.ChargeCarrier(
        name="e", charge=-1, mobility=mobility, diffusivity=0.04
    )
    positive_ions = ic.ChargeCarrier(name="P", charge=1, mobility=2.0e-4)
    recombination = ic.CarrierReaction(
        reactant=[positive_ions, electrons], rate=2.0e-13
    )
    return ic.Gas(
        name="non-attaching",
        w_value=33.97,
        density=1.205,
        carriers=[electrons, positive_ions],
        ionisation_products=(electrons, positive_ions),
        reactions=[recombination],
    )


def run_coarse_problem(quantity, gas=None, bias=500.0, activity=3.7e4) -> float:
    """Solves the reference chamber on a coarse mesh without wall losses and returns
    the value of one quantity."""
    chamber = reference_chamber()
    if isinstance(quantity, str):
        # the quantity needs the anode of this chamber
        quantity = getattr(ic, quantity)(surface=chamber.anode)
    transport = ic.ChargeTransportProblem(
        electrostatic_problem=solve_potential(
            chamber, radial_cells=30, axial_cells=18, bias=bias
        ),
        gas=gas or air(),
        ionisation_sources=[ic.TritiumInGas(activity=activity)],
        exports=[quantity],
    )
    transport.initialise()
    transport.run()
    return quantity.value


# ---- structure ----------------------------------------------------------------------


def test_is_a_festim_hydrogen_transport_problem():
    """A ChargeTransportProblem is a FESTIM HydrogenTransportProblem."""
    assert isinstance(ic.ChargeTransportProblem(), F.HydrogenTransportProblem)


def test_chamber_quantities_are_festim_derived_quantities():
    """A chamber quantity is a FESTIM derived quantity."""
    assert isinstance(ic.CollectionEfficiency(), F.exports.DerivedQuantity)


def test_initialise_needs_a_solved_electrostatic_problem(chamber):
    """Initialising with an electrostatic problem that was never run raises a
    RuntimeError."""
    transport = ic.ChargeTransportProblem(
        electrostatic_problem=ic.ElectrostaticProblem(), gas=air()
    )
    with pytest.raises(RuntimeError, match="electrostatic problem must be"):
        transport.initialise()


def test_uses_the_mesh_of_the_electrostatic_problem(coarse_transport):
    """The transport problem is solved on the mesh of the electrostatic problem."""
    coarse_transport.initialise()
    assert coarse_transport.mesh is coarse_transport.electrostatic_problem.mesh


def test_volume_of_gas_is_the_chamber_volume(coarse_transport, chamber):
    """The gas volume integrated over the mesh equals the volume of the chamber."""
    coarse_transport.initialise()
    assert np.isclose(coarse_transport.volume_of_gas, chamber.volume)


def test_reduced_field_of_the_reference_chamber(coarse_transport):
    """The largest E/N in the reference chamber at 500 V lies between 1 and 10 Td."""
    coarse_transport.initialise()
    assert 1.0 < coarse_transport.reduced_field.max() < 10.0


# ---- stabilisation ------------------------------------------------------------------


@pytest.mark.parametrize(
    "diffusivity, expected",
    [(1e-6, 5e-4), (1.0, 1.0)],  # floor active, floor inactive
)
def test_floor_is_the_larger_of_physical_and_upwind(diffusivity, expected):
    """The floored diffusivity is the upwind value when that is larger and the
    physical value otherwise."""
    D = ic.DiffusionFloor().apply(
        diffusivity=np.array([diffusivity]),
        mobility=np.array([1e-4]),
        field_strength=np.array([1e4]),
        cell_length=np.array([1e-3]),
    )
    assert np.isclose(D[0], expected)


def test_floor_factor_zero_switches_it_off():
    """A floor factor of zero leaves the physical diffusivity unchanged."""
    D = ic.DiffusionFloor(factor=0.0).apply(
        np.array([1e-6]), np.array([1e-4]), np.array([1e4]), np.array([1e-3])
    )
    assert D[0] == 1e-6


def test_floor_rejects_negative_factor():
    """A negative floor factor raises a ValueError."""
    with pytest.raises(ValueError, match="factor must not be negative"):
        ic.DiffusionFloor(factor=-1.0)


# ---- regression against the prototype (chamber_model.py, N_r = 60) ------------------


@pytest.mark.parametrize(
    "index, expected",
    [
        (MICROCURIE, 9.40662259521735e-16),
        (CURIE, 9.33952139344237e-10),
        (THOUSAND_CURIE, 3.2356679629878293e-07),
    ],
)
def test_current_matches_the_prototype(reference_results, index, expected):
    """The current of the reference chamber equals the original prototype's to 8
    digits."""
    assert np.isclose(reference_results["current"].data[index], expected, rtol=1e-8)


def test_wall_efficiency_matches_the_prototype(reference_results):
    """The wall efficiency of the reference chamber equals the original prototype's
    to 8 digits."""
    value = reference_results["wall_efficiency"].data[MICROCURIE]
    assert np.isclose(value, 0.9895863173320661, rtol=1e-8)


def test_collection_efficiency_matches_the_prototype(reference_results):
    """The collection efficiency at 1000 Ci/m3 equals the original prototype's to 8
    digits."""
    value = reference_results["collection_efficiency"].data[THOUSAND_CURIE]
    assert np.isclose(value, 0.36235411253601335, rtol=1e-8)


def test_space_charge_field_matches_the_prototype(reference_results):
    """The space-charge field ratio at 1000 Ci/m3 equals the original prototype's to
    8 digits."""
    value = reference_results["space_charge_field"].data[THOUSAND_CURIE]
    assert np.isclose(value, 0.6501553934189711, rtol=1e-8)


# ---- physics checks -----------------------------------------------------------------


@pytest.mark.parametrize("index", [MICROCURIE, CURIE, THOUSAND_CURIE])
def test_ramo_current_agrees_with_the_anode_flux(reference_results, index):
    """The Shockley-Ramo current and the charge collected by the anode agree within
    0.1%. The anode collects negative charge, hence the minus sign."""
    ramo = reference_results["current"].data[index]
    anode = reference_results["anode_current"].data[index]
    assert np.isclose(ramo, -anode, rtol=1e-3)


@pytest.mark.parametrize("index", [MICROCURIE, CURIE])
def test_charge_leaving_through_both_electrodes_cancels(reference_results, index):
    """The charge leaving through the anode and through the cathode are opposite
    within 0.5%, as charge conservation requires."""
    anode = reference_results["anode_current"].data[index]
    cathode = reference_results["cathode_current"].data[index]
    assert np.isclose(anode, -cathode, rtol=5e-3)


def test_response_is_linear_at_low_activity(reference_results):
    """Ten times the activity gives ten times the current while recombination is
    negligible."""
    currents = reference_results["current"].data
    assert np.isclose(currents[TEN_MICROCURIE], 10 * currents[MICROCURIE], rtol=1e-5)


def test_recombination_is_negligible_at_low_activity(reference_results):
    """At 1 uCi/m3 fewer than one ion pair in a million recombines."""
    assert reference_results["collection_efficiency"].data[MICROCURIE] > 1 - 1e-6


def test_recombination_dominates_at_high_activity(reference_results):
    """At 1000 Ci/m3 more than half of the ion pairs recombine."""
    assert reference_results["collection_efficiency"].data[THOUSAND_CURIE] < 0.5


def test_space_charge_grows_with_activity(reference_results):
    """The space-charge field ratio increases with the activity."""
    fields = reference_results["space_charge_field"].data
    assert fields == sorted(fields)


def test_current_never_exceeds_the_ideal_current(reference_results):
    """The current stays below the ideal current at every activity."""
    currents = np.array(reference_results["current"].data)
    ideal = np.array(reference_results["ideal_current"].data)
    assert np.all(currents < ideal)


def test_collection_efficiency_is_consistent_with_the_rates(reference_results):
    """The collection efficiency equals (generation - recombination) / generation."""
    generation = reference_results["generation"].data[CURIE]
    recombination = reference_results["recombination"].data[CURIE]
    expected = (generation - recombination) / generation
    value = reference_results["collection_efficiency"].data[CURIE]
    assert np.isclose(value, expected)


def test_data_holds_one_value_per_run(reference_results):
    """A quantity keeps the value of every run in data."""
    assert len(reference_results["current"].data) == len(ACTIVITIES)


def test_without_wall_losses_every_pair_is_generated():
    """A source without a geometry gives a wall efficiency of one."""
    assert np.isclose(run_coarse_problem(ic.WallEfficiency()), 1.0)


def test_reversing_the_bias_reverses_the_current():
    """A negative bias gives the same current with the opposite sign."""
    forward = run_coarse_problem("ChamberCurrent", bias=500.0)
    reverse = run_coarse_problem("ChamberCurrent", bias=-500.0)
    assert np.isclose(reverse, -forward, rtol=1e-6)


def test_higher_bias_collects_more():
    """Raising the bias from 100 V to 1000 V raises the collection efficiency."""
    low = run_coarse_problem(ic.CollectionEfficiency(), bias=100.0, activity=3.7e12)
    high = run_coarse_problem(ic.CollectionEfficiency(), bias=1000.0, activity=3.7e12)
    assert high > low


def test_thinner_gas_loses_more_to_the_walls(chamber):
    """Air at a fifth of its normal density has a lower wall efficiency."""
    efficiencies = []
    for density in (1.205, 0.241):
        gas = air()
        gas.density = density
        wall_efficiency = ic.WallEfficiency()
        transport = ic.ChargeTransportProblem(
            electrostatic_problem=solve_potential(
                chamber, radial_cells=30, axial_cells=18
            ),
            gas=gas,
            ionisation_sources=[
                ic.TritiumInGas(activity=3.7e4, geometry=chamber, n_rays=300)
            ],
            exports=[wall_efficiency],
        )
        transport.initialise()
        transport.run()
        efficiencies.append(wall_efficiency.value)
    assert efficiencies[1] < efficiencies[0]


def test_chamber_current_needs_an_electrode_of_the_problem(coarse_transport):
    """A ChamberCurrent on a surface that is not an electrode raises a ValueError."""
    stranger = F.SurfaceSubdomain(id=9)
    coarse_transport.exports = [ic.ChamberCurrent(surface=stranger)]
    coarse_transport.initialise()
    with pytest.raises(ValueError, match="is not an electrode"):
        coarse_transport.run()


def test_quantity_writes_a_csv_file(coarse_transport, tmp_path):
    """A quantity given a filename writes its value to that file."""
    filename = tmp_path / "collection.csv"
    coarse_transport.exports = [ic.CollectionEfficiency(filename=str(filename))]
    coarse_transport.initialise()
    coarse_transport.run()
    assert filename.exists()


def test_festim_species_export_writes_a_file(coarse_transport, tmp_path):
    """FESTIM's own SpeciesExport works on a transport problem and writes the carrier
    density to a file."""
    filename = tmp_path / "electrons.bp"
    electrons = coarse_transport.gas.carriers[0]
    coarse_transport.exports = [F.SpeciesExport(field=electrons, filename=filename)]
    coarse_transport.initialise()
    coarse_transport.run()
    assert filename.exists()


# ---- flexibility of the gas description ---------------------------------------------


def test_gas_without_negative_ions_collects_nearly_everything():
    """A gas of electrons and positive ions only solves and collects within 10% of
    the ideal current."""
    chamber = reference_chamber()
    current = ic.ChamberCurrent(surface=chamber.anode)
    ideal = ic.IdealCurrent()
    transport = ic.ChargeTransportProblem(
        electrostatic_problem=solve_potential(chamber, radial_cells=30, axial_cells=18),
        gas=electron_ion_gas(),
        ionisation_sources=[ic.TritiumInGas(activity=3.7e4)],
        exports=[current, ideal],
    )
    transport.initialise()
    transport.run()
    assert np.isclose(current.value, ideal.value, rtol=0.1)


def test_gas_without_reactions_has_no_recombination():
    """A gas with no reactions reports zero recombination."""
    gas = electron_ion_gas()
    gas.reactions = []
    assert run_coarse_problem(ic.Recombination(), gas=gas) == 0.0


def test_callable_mobility_matches_the_equivalent_constant():
    """A mobility function returning a constant gives the same current as that
    constant."""

    def constant_mobility(reduced_field):
        return np.full_like(reduced_field, 0.04)

    with_float = run_coarse_problem("ChamberCurrent", gas=electron_ion_gas(0.04))
    with_callable = run_coarse_problem(
        "ChamberCurrent", gas=electron_ion_gas(constant_mobility)
    )
    assert np.isclose(with_float, with_callable, rtol=1e-9)


def test_two_sources_add_up(coarse_electrostatics):
    """Two sources generate as many ion pairs as one source with their summed
    activity."""
    generation = ic.PairGeneration()
    ideal = ic.IdealCurrent()
    transport = ic.ChargeTransportProblem(
        electrostatic_problem=coarse_electrostatics,
        gas=air(),
        ionisation_sources=[
            ic.TritiumInGas(activity=1e4),
            ic.TritiumInGas(activity=3e4),
        ],
        exports=[generation, ideal],
    )
    transport.initialise()
    transport.run()
    one_source = ic.TritiumInGas(activity=4e4).ideal_pair_rate(
        air(), transport.volume_of_gas
    )
    assert np.isclose(generation.value, one_source)


# ---- a mesh read from a file --------------------------------------------------------

GAS_ID = 7
ANODE_ID = 11
CATHODE_ID = 12


def write_tagged_mesh(chamber: ic.CoaxialChamber, filename) -> None:
    """Writes the mesh of a chamber to a file with cell and facet tags, as a mesh
    generator would."""
    mesh = chamber.create_mesh(radial_cells=30, axial_cells=18).mesh

    num_cells = mesh.topology.index_map(2).size_local
    all_cells = np.arange(num_cells, dtype=np.int32)
    cell_values = np.full(num_cells, GAS_ID, dtype=np.int32)
    cell_tags = dolfinx.mesh.meshtags(mesh, 2, all_cells, cell_values)

    anode_facets = dolfinx.mesh.locate_entities_boundary(mesh, 1, chamber.is_anode)
    cathode_facets = dolfinx.mesh.locate_entities_boundary(mesh, 1, chamber.is_cathode)
    facets = np.concatenate([anode_facets, cathode_facets])
    facet_values = np.concatenate(
        [
            np.full(len(anode_facets), ANODE_ID, dtype=np.int32),
            np.full(len(cathode_facets), CATHODE_ID, dtype=np.int32),
        ]
    )
    order = np.argsort(facets)
    facet_tags = dolfinx.mesh.meshtags(mesh, 1, facets[order], facet_values[order])

    io4dolfinx.write_mesh(filename, mesh)
    io4dolfinx.write_meshtags(filename, mesh, cell_tags, meshtag_name="cells")
    io4dolfinx.write_meshtags(filename, mesh, facet_tags, meshtag_name="facets")


def run_on_mesh_from_file(filename) -> float:
    """Reads a mesh and its tags from a file, solves both problems on it with
    subdomains defined by their ids only, and returns the chamber current."""
    mesh = io4dolfinx.read_mesh(filename, MPI.COMM_WORLD)
    cell_tags = io4dolfinx.read_meshtags(filename, mesh, meshtag_name="cells")
    facet_tags = io4dolfinx.read_meshtags(filename, mesh, meshtag_name="facets")

    gas_volume = F.VolumeSubdomain(id=GAS_ID, material=None)
    anode = F.SurfaceSubdomain(id=ANODE_ID)
    cathode = F.SurfaceSubdomain(id=CATHODE_ID)

    electrostatics = ic.ElectrostaticProblem()
    electrostatics.mesh = F.Mesh(mesh, coordinate_system="cylindrical")
    electrostatics.volume_meshtags = cell_tags
    electrostatics.facet_meshtags = facet_tags
    electrostatics.subdomains = [gas_volume, anode, cathode]
    electrostatics.boundary_conditions = [
        ic.FixedPotentialBC(subdomain=anode, value=500.0),
        ic.FixedPotentialBC(subdomain=cathode, value=0.0),
    ]
    electrostatics.initialise()
    electrostatics.run()

    current = ic.ChamberCurrent(surface=anode)
    transport = ic.ChargeTransportProblem(
        electrostatic_problem=electrostatics,
        gas=air(),
        ionisation_sources=[ic.TritiumInGas(activity=3.7e4)],
        exports=[current, ic.SpaceChargeField()],
    )
    transport.initialise()
    transport.run()
    return current.value


def test_mesh_read_from_a_file_gives_the_same_current(tmp_path):
    """A mesh read from a file, with subdomains defined by the ids of its tags, gives
    the same current as the same mesh created by the chamber."""
    filename = tmp_path / "chamber.bp"
    write_tagged_mesh(reference_chamber(), filename)
    expected = run_coarse_problem("ChamberCurrent")
    assert np.isclose(run_on_mesh_from_file(filename), expected, rtol=1e-8)


def test_gas_volume_takes_the_id_of_the_electrostatic_problem(coarse_transport):
    """After initialise, the gas volume of the transport problem has the id of the
    volume subdomain of the electrostatic problem."""
    coarse_transport.initialise()
    expected = coarse_transport.electrostatic_problem.volume_subdomains[0].id
    assert coarse_transport.gas_volume.id == expected


def test_rejects_an_electrostatic_problem_with_two_volumes(chamber):
    """An electrostatic problem with more than one volume subdomain raises a
    ValueError, because the transport problem solves in the gas only."""
    electrostatics = solve_potential(chamber, radial_cells=30, axial_cells=18)
    electrostatics.subdomains.append(F.VolumeSubdomain(id=2, material=None))
    transport = ic.ChargeTransportProblem(
        electrostatic_problem=electrostatics,
        gas=air(),
        ionisation_sources=[ic.TritiumInGas(activity=3.7e4)],
    )
    with pytest.raises(ValueError, match="must have one volume subdomain"):
        transport.initialise()


def test_rejects_a_mesh_that_is_not_made_of_triangles():
    """A mesh of quadrilaterals raises a ValueError, because the stabilisation is
    written for triangles."""
    quadrilaterals = dolfinx.mesh.create_unit_square(
        MPI.COMM_WORLD, 4, 4, dolfinx.mesh.CellType.quadrilateral
    )
    left = F.SurfaceSubdomain(id=1, locator=lambda x: np.isclose(x[0], 0.0))
    right = F.SurfaceSubdomain(id=2, locator=lambda x: np.isclose(x[0], 1.0))
    electrostatics = ic.ElectrostaticProblem(
        mesh=F.Mesh(quadrilaterals),
        subdomains=[F.VolumeSubdomain(id=1, material=None), left, right],
        boundary_conditions=[
            ic.FixedPotentialBC(subdomain=left, value=100.0),
            ic.FixedPotentialBC(subdomain=right, value=0.0),
        ],
    )
    electrostatics.initialise()
    electrostatics.run()
    transport = ic.ChargeTransportProblem(
        electrostatic_problem=electrostatics,
        gas=air(),
        ionisation_sources=[ic.TritiumInGas(activity=3.7e4)],
    )
    with pytest.raises(ValueError, match="must be made of triangles"):
        transport.initialise()
