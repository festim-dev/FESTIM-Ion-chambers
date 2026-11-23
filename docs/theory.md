# Theory

## Charge transport

Each charge carrier $s$ has a number density $c_s$ (m⁻³) that obeys a steady-state
Nernst–Planck balance:

$$
\nabla \cdot J_s = S_s + R_s, \qquad
J_s = -D_s \nabla c_s + c_s v_s, \qquad
v_s = -z_s\, \mu_s \nabla \varphi
$$

with $D_s$ the diffusion coefficient, $\mu_s$ the mobility, $z_s$ the charge number
and $\varphi$ the electric potential. $S_s$ is the ionisation source and $R_s$ the net
production by reactions. The drift term is FESTIM's `ElectromigrationTerm`, with the
velocity set by the mobility
({class}`~ion_chambers.MobilityElectromigration`). This is what
{class}`~ion_chambers.ChargeTransportProblem` solves.

For air, the three balances are

$$
\begin{aligned}
\nabla \cdot J_e &= S - \nu_a n_e - \beta\, n_P n_e \\
\nabla \cdot J_P &= S - \alpha\, n_P n_N - \beta\, n_P n_e \\
\nabla \cdot J_N &= \nu_a n_e - \alpha\, n_P n_N
\end{aligned}
$$

for electrons, positive ions and negative ions, with $\nu_a$ the attachment frequency,
$\alpha$ the ion–ion and $\beta$ the electron–ion recombination coefficient. In air,
electrons attach to O₂ within about 10 ns, a few micrometres from where they are
created, so nearly all the negative charge travels as slow ions.

## Electric potential

The potential is the solution of an {class}`~ion_chambers.ElectrostaticProblem`. It
solves Laplace's equation, $\nabla^2 \varphi = 0$, with $\varphi = V$ on
the anode, $\varphi = 0$ on the cathode and no normal field on the insulators. The
charge in the gas is not fed back into the potential.

To check that assumption, the potential $\varphi_1$ created by the converged charge
density is computed afterwards from

$$
-\nabla \cdot (\varepsilon_0 \nabla \varphi_1) = e \sum_s z_s c_s,
\qquad \varphi_1 = 0 \text{ on the electrodes}
$$

and the ratio $\lVert \nabla\varphi_1 \rVert / \lVert \nabla\varphi \rVert$ is
reported by the export {class}`~ion_chambers.SpaceChargeField`.

## Source

Tritium mixed into the gas at activity concentration $A$ (Bq/m³) creates ion pairs at
the rate

$$
S(r, z) = A\, \frac{\bar{E}}{W}\, f(r, z)
$$

where $\bar{E}$ = 5.70 keV is the mean beta energy, $W$ the mean energy per ion pair in
the gas (33.97 eV in air, giving 167.7 pairs per decay) and $f$ the wall-loss
efficiency.

### Wall losses

$f(r, z)$ is the share of the beta energy that a decay at $(r, z)$ deposits in the
gas. It is computed by a Monte Carlo with these assumptions:

- Energies follow the tritium beta spectrum (allowed shape with the Fermi function,
  endpoint 18.59 keV).
- Emission is isotropic.
- Tracks are straight, with the Katz–Penfold practical range: 4.9 mm in air at the
  endpoint energy. Scattering is folded into the practical range.
- A beta that reaches a wall deposits there the energy it still has, following the
  range–energy relation.
- The energy is counted at the point of decay. The track is short compared with the
  drift gap.

## Boundary conditions

An electrode absorbs the carriers the field drives into it (an outflow condition) and
emits none (zero density for the carriers the field pushes away). Insulators carry no
flux.

## Chamber current

The reading is the current induced on the anode, from the Shockley–Ramo theorem:

$$
I = e \int \sum_s z_s\, J_s \cdot \left(-\nabla \varphi_w\right) dV,
\qquad \varphi_w = \varphi / V
$$

where $\varphi_w$ is the weighting potential. In steady state this equals the current
through either electrode.

## What the model leaves out

- **Space-charge feedback.** The field is fixed. The distortion is reported but not
  applied.
- **Gas flow.** The gas is still. In a flow-through chamber the flow is two to four
  orders of magnitude slower than the ion drift.
- **Columnar recombination** along individual beta tracks. Only volume recombination
  between uniformly mixed carriers is modelled.
- **Tritium on the walls**, which causes memory effects and background drift.
- **Time dependence.** Only the steady state is solved.
