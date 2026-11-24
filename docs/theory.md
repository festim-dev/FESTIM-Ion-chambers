# Theory

An ionisation chamber turns radioactive decays in a gas into an electric current. This
page follows what happens between the two, one mechanism at a time, and gives the
equation the model uses for each. The [walkthrough](walkthrough.md) builds the same
model in code.

| | Mechanism | What happens | Scale |
|---|---|---|---|
| 1 | [Beta decay](#beta-decay) | A tritium nucleus emits a fast electron, the beta. | 5.7 keV on average |
| 2 | [Ionisation](#ionisation) | The beta knocks electrons out of gas molecules: ion pairs. | 168 pairs within about 0.4 mm |
| 3 | [Wall losses](#wall-losses) | A beta emitted near a wall ends in the metal. | within 5 mm of a wall |
| 4 | [Attachment](#attachment) | A free electron sticks to O₂ and becomes a negative ion. | after 10 ns |
| 5 | [Drift and diffusion](#drift-and-diffusion) | The field pulls positive and negative charges apart. | 0.6 to 19 m/s |
| 6 | [Recombination](#recombination) | A positive and a negative ion meet and neutralise. | grows with the activity |
| 7 | [Collection](#collection-and-the-chamber-current) | The charges reach the electrodes: the current. | 0.97 nA at 1 Ci/m³ |

The numbers on this page are for the chamber of the [walkthrough](walkthrough.md)
(very nearly 1 L) in air at 500 V, with the
illustrative gas data of the [example case](example_case.md#the-script).

## Beta decay

```{raw} html
<div class="ic-figure">
<svg viewBox="0 0 470 190" role="img" aria-label="A tritium nucleus turns into helium-3 and emits a beta and an antineutrino">
        <circle cx="62" cy="86" r="13" fill="var(--ic-pos)"/>
        <circle cx="84" cy="98" r="13" fill="var(--ic-mol)"/>
        <circle cx="66" cy="110" r="13" fill="var(--ic-mol)"/>
        <text x="72" y="150" font-size="13" fill="var(--ic-ink)" text-anchor="middle">tritium</text>
        <text x="72" y="167" font-size="12" fill="var(--ic-muted)" text-anchor="middle">1 proton, 2 neutrons</text>
        <path d="M122 98 L186 98 M178 92 L186 98 L178 104" stroke="var(--ic-ink)" stroke-width="2" fill="none"/>
        <circle cx="232" cy="86" r="13" fill="var(--ic-pos)"/>
        <circle cx="254" cy="98" r="13" fill="var(--ic-pos)"/>
        <circle cx="236" cy="110" r="13" fill="var(--ic-mol)"/>
        <text x="242" y="150" font-size="13" fill="var(--ic-ink)" text-anchor="middle">helium-3</text>
        <text x="242" y="167" font-size="12" fill="var(--ic-muted)" text-anchor="middle">2 protons, 1 neutron</text>
        <path d="M276 84 L350 44 M340 44 L350 44 L346 53" stroke="var(--ic-beta)" stroke-width="2" fill="none"/>
        <circle cx="362" cy="38" r="6" fill="var(--ic-beta)"/>
        <text x="374" y="30" font-size="13" fill="var(--ic-ink)">beta</text>
        <text x="374" y="46" font-size="12" fill="var(--ic-muted)">an electron</text>
        <path d="M278 108 L352 136" stroke="var(--ic-muted)" stroke-width="1.5" stroke-dasharray="3 4" fill="none"/>
        <text x="358" y="134" font-size="12" fill="var(--ic-muted)">antineutrino</text>
        <text x="358" y="150" font-size="12" fill="var(--ic-muted)">escapes unseen</text>
      </svg>
</div>
```

One of the two neutrons of the tritium nucleus turns into a proton. The nucleus is now
helium-3, and it throws out a fast electron, the beta, together with an antineutrino
that leaves the chamber without touching anything.

- The decay releases 18.6 keV. The beta and the antineutrino share it at random, so a
  beta has anything from 0 to 18.6 keV. The average is $\bar{E}$ = 5.70 keV.
- Decays happen at random. Half of any amount of tritium decays in 12.3 years.
- The activity counts them: 1 Bq is one decay per second, and 1 Ci is 3.7 × 10¹⁰ of
  them.

**In the model.** The energies follow the tritium beta spectrum: the allowed shape
with the Fermi function and an endpoint of 18.59 keV
({func}`ion_chambers.beta.beta_spectrum`).

## Ionisation

The beta gives a little of its energy to an electron of a gas molecule it passes. If
that is enough to pull the electron off, the molecule is left as a positive ion and the
electron is free. Those two are one *ion pair*. The beta carries on, a little slower,
and does it again.

```{raw} html
<div class="ic-steps">
<div class="ic-step">
<svg viewBox="0 0 280 180" role="img" aria-label="A beta approaches a neutral gas molecule">
          <ellipse cx="150" cy="104" rx="44" ry="30" fill="var(--ic-mol)" opacity="0.35"/>
          <circle cx="134" cy="104" r="9" fill="var(--ic-mol)"/>
          <circle cx="166" cy="104" r="9" fill="var(--ic-mol)"/>
          <circle cx="150" cy="78" r="4" fill="var(--ic-neg)"/>
          <circle cx="118" cy="118" r="4" fill="var(--ic-neg)"/>
          <circle cx="184" cy="120" r="4" fill="var(--ic-neg)"/>
          <line x1="20" y1="52" x2="84" y2="52" stroke="var(--ic-beta)" stroke-width="2" stroke-dasharray="4 4"/>
          <circle cx="92" cy="52" r="6" fill="var(--ic-beta)"/>
          <path d="M104 52 L122 52 M116 47 L122 52 L116 57" stroke="var(--ic-beta)" stroke-width="2" fill="none"/>
          <text x="20" y="34" font-size="13" fill="var(--ic-ink)">beta, 5 700 eV</text>
          <text x="150" y="160" font-size="13" fill="var(--ic-muted)" text-anchor="middle">neutral molecule (N₂)</text>
        </svg>
<strong>1. The beta comes close</strong>
The beta is a fast electron. A gas molecule is two nuclei surrounded by its own electrons, and has no net charge.
</div>
<div class="ic-step">
<svg viewBox="0 0 280 180" role="img" aria-label="The beta pushes one of the molecule's electrons away">
          <ellipse cx="150" cy="104" rx="44" ry="30" fill="var(--ic-mol)" opacity="0.35"/>
          <circle cx="134" cy="104" r="9" fill="var(--ic-mol)"/>
          <circle cx="166" cy="104" r="9" fill="var(--ic-mol)"/>
          <circle cx="118" cy="118" r="4" fill="var(--ic-neg)"/>
          <circle cx="184" cy="120" r="4" fill="var(--ic-neg)"/>
          <circle cx="150" cy="52" r="6" fill="var(--ic-beta)"/>
          <line x1="78" y1="52" x2="140" y2="52" stroke="var(--ic-beta)" stroke-width="2" stroke-dasharray="4 4"/>
          <circle cx="156" cy="84" r="4" fill="var(--ic-neg)"/>
          <path d="M152 60 L155 74" stroke="var(--ic-ink)" stroke-width="1.5" stroke-dasharray="2 3" fill="none"/>
          <path d="M160 92 L176 124 M168 120 L176 124 L177 115" stroke="var(--ic-neg)" stroke-width="2" fill="none"/>
          <text x="166" y="48" font-size="13" fill="var(--ic-ink)">like charges repel</text>
          <text x="150" y="160" font-size="13" fill="var(--ic-muted)" text-anchor="middle">one electron gets a kick</text>
        </svg>
<strong>2. It pushes on one electron</strong>
Both are negative, so they repel. The beta is gone in an instant, but the push leaves that electron with extra energy.
</div>
<div class="ic-step">
<svg viewBox="0 0 280 180" role="img" aria-label="A positive ion and a free electron are left behind">
          <ellipse cx="120" cy="104" rx="44" ry="30" fill="var(--ic-pos)" opacity="0.25"/>
          <circle cx="104" cy="104" r="9" fill="var(--ic-pos)"/>
          <circle cx="136" cy="104" r="9" fill="var(--ic-pos)"/>
          <circle cx="88" cy="118" r="4" fill="var(--ic-neg)"/>
          <circle cx="154" cy="120" r="4" fill="var(--ic-neg)"/>
          <text x="120" y="88" font-size="18" font-weight="600" fill="var(--ic-pos)" text-anchor="middle">+</text>
          <circle cx="218" cy="132" r="4" fill="var(--ic-neg)"/>
          <path d="M184 116 L208 128" stroke="var(--ic-neg)" stroke-width="1.5" stroke-dasharray="3 3" fill="none"/>
          <text x="228" y="150" font-size="13" fill="var(--ic-neg)" text-anchor="middle">free electron</text>
          <line x1="150" y1="50" x2="226" y2="44" stroke="var(--ic-beta)" stroke-width="2" stroke-dasharray="4 4"/>
          <circle cx="234" cy="43" r="6" fill="var(--ic-beta)"/>
          <text x="150" y="30" font-size="13" fill="var(--ic-ink)">beta, about 5 665 eV</text>
          <text x="112" y="160" font-size="13" fill="var(--ic-pos)" text-anchor="middle">positive ion (N₂⁺)</text>
        </svg>
<strong>3. An ion pair is left</strong>
The electron leaves. The molecule, one electron short, is now a positive ion. The beta has lost what it gave and moves on.
</div>
</div>
```

The beta never hits a nucleus and nothing is destroyed. It only moves an electron from
a molecule to the open gas. Before, the gas could carry no current. Now it holds two
charges that the electric field of the chamber can pull apart.

The beta itself is far too fast for the field of the chamber to steer. Only its
collisions with the gas deflect it, until its energy is spent. Then it is one more slow electron
among those it freed.

### The track of one beta

Follow one beta of 5 700 eV until it stops. Not every molecule it passes is ionised.

```{raw} html
<div class="ic-stage">
<div class="ic-readouts">
<div class="ic-readout"><span class="ic-label">Energy left</span><span class="ic-value" id="ic-energy">5 700 eV</span><div class="ic-meter"><span id="ic-meter"></span></div></div>
<div class="ic-readout"><span class="ic-label">Ion pairs</span><span class="ic-value" id="ic-pairs">0</span></div>
<div class="ic-readout"><span class="ic-label">Spent per pair</span><span class="ic-value" id="ic-per-pair">–</span></div>
</div>
<canvas id="ic-track" aria-label="Animation of a beta crossing the gas and leaving ion pairs behind"></canvas>
<div class="ic-controls">
<button id="ic-play" type="button">Pause</button>
<button id="ic-restart" type="button">New beta</button>
<div class="ic-group" role="group" aria-label="Speed"><span class="ic-label">Speed</span>
<button id="ic-speed1" type="button" aria-pressed="false">1×</button>
<button id="ic-speed4" type="button" aria-pressed="true">4×</button>
<button id="ic-speed16" type="button" aria-pressed="false">16×</button></div>
<span class="ic-status" id="ic-status" role="status"></span>
</div>
</div>
<div class="ic-legend">
<span><i class="ic-dot" style="background: var(--ic-beta)"></i>the beta</span>
<span><i class="ic-dot" style="background: var(--ic-mol)"></i>neutral molecule</span>
<span><i class="ic-dot" style="background: var(--ic-pos)"></i>positive ion</span>
<span><i class="ic-dot" style="background: var(--ic-neg)"></i>free electron</span>
<span><i class="ic-ring"></i>excited molecule</span>
</div>
```

:::{note}
The animation is not to scale. Real molecules are about 0.3 nm across and 3.4 nm
apart, and a beta of the average energy ends about 0.39 mm from where it started. It leaves
one pair every few micrometres, so it passes hundreds of molecules untouched for each
one it ionises. The
energies are drawn at random so that the average matches air.
:::

Three things can happen at a molecule:

| Outcome | What the beta gives | Result |
|---|---|---|
| Nothing | 0 eV | The beta passes too far away. This is by far the most common. |
| Excitation | about 10 eV | The push is too weak to free an electron. The molecule gives the energy back as light or heat. No charge is made. |
| Ionisation | 15.6 eV or more | An electron is freed: one ion pair. Anything above 15.6 eV becomes the speed of that electron. |

Sometimes the freed electron is itself fast enough to ionise another molecule nearby.
These count as ion pairs too.

### The energy per ion pair

Pulling an electron off a nitrogen molecule takes 15.6 eV. Yet the beta spends 34 eV
on average for each pair it leaves in air. The rest goes into excitations and into the
speed of the freed electrons, and makes no charge.

This average is called $W$. It belongs to the gas: about 34 eV in air, 35 eV in
nitrogen and 26 eV in argon. It barely depends on the energy of the beta, so a decay
makes $\bar{E} / W$ pairs: 168 in air.

**In the model.** The collisions are not followed. Tritium mixed into the gas at the
activity concentration $A$ (Bq/m³) creates ion pairs at the rate

$$
S(r, z) = A\, \frac{\bar{E}}{W}\, f(r, z)
$$

where $f$ is the wall-loss efficiency of the next section
({class}`~ion_chambers.TritiumInGas`). The pairs are created at the point of the
decay, because the track (0.39 mm at the average energy) is short next to the gas gap
(48.5 mm).

## Wall losses

A beta goes a set distance before it stops, its *range*. If a wall is closer than
that, the beta ends in the metal, and the energy it still had makes no ion pairs. Move
the decay towards the wall to see how much is lost.

```{raw} html
<div class="ic-stage">
<div class="ic-readouts">
<div class="ic-readout"><span class="ic-label">Betas that reach the wall</span><span class="ic-value" id="ic-wall-hits">–</span></div>
<div class="ic-readout"><span class="ic-label">Energy left in the gas</span><span class="ic-value" id="ic-wall-kept">–</span></div>
<div class="ic-readout"><span class="ic-label">Longest possible track</span><span class="ic-value" id="ic-wall-range">–</span></div>
</div>
<canvas id="ic-wall-canvas" aria-label="Tracks of betas from a decay near a wall; those that reach the wall are red"></canvas>
<div class="ic-controls">
<label class="ic-slider" for="ic-wall-distance">Distance from the wall <output id="ic-wall-distance-out">1.00 mm</output>
<input id="ic-wall-distance" type="range" min="0" max="5" step="0.05" value="1"></label>
<label class="ic-slider" for="ic-wall-pressure">Gas pressure (air) <output id="ic-wall-pressure-out">1.0 atm</output>
<input id="ic-wall-pressure" type="range" min="0.2" max="2" step="0.1" value="1"></label>
</div>
</div>
<div class="ic-legend">
<span><i class="ic-dot" style="background: var(--ic-muted)"></i>stops in the gas</span>
<span><i class="ic-dot" style="background: var(--ic-pos)"></i>reaches the wall</span>
</div>
```

- **The range depends on the energy.** A beta of the average energy goes 0.39 mm in
  air. The rare ones near 18.6 keV go 4.9 mm. So only decays within a few millimetres
  of a wall lose anything.
- **The range depends on the gas density.** A beta stops after crossing a fixed mass
  of gas. Halve the pressure and every track is twice as long. A light gas such as
  helium does the same.
- **Any wall counts**: the can, the end caps and the rod.

**In the model.** $f(r, z)$ is the share of the beta energy that a decay at $(r, z)$
deposits in the gas. It is computed by a Monte Carlo
({func}`ion_chambers.beta.wall_efficiency`) with these assumptions:

- Energies follow the tritium beta spectrum.
- Emission is isotropic.
- Tracks are straight, with the Katz–Penfold practical range: 4.9 mm in air at the
  endpoint energy. Scattering is folded into the practical range.
- A beta that reaches a wall deposits there the energy it still has, following the
  range–energy relation.
- The energy is counted at the point of decay.

The drawing above repeats this calculation for a flat wall. In the reference chamber,
1.1% of the beta energy is lost to the walls.

## Attachment

```{raw} html
<div class="ic-figure">
<svg viewBox="0 0 440 120" role="img" aria-label="A free electron attaches to an oxygen molecule, which becomes a negative ion">
        <circle cx="50" cy="56" r="5" fill="var(--ic-neg)"/>
        <text x="50" y="96" font-size="12" fill="var(--ic-muted)" text-anchor="middle">free electron</text>
        <text x="96" y="62" font-size="18" fill="var(--ic-ink)" text-anchor="middle">+</text>
        <circle cx="140" cy="56" r="13" fill="var(--ic-mol)"/>
        <circle cx="162" cy="56" r="13" fill="var(--ic-mol)"/>
        <text x="151" y="96" font-size="12" fill="var(--ic-muted)" text-anchor="middle">oxygen molecule</text>
        <path d="M204 56 L262 56 M254 50 L262 56 L254 62" stroke="var(--ic-ink)" stroke-width="2" fill="none"/>
        <circle cx="312" cy="56" r="13" fill="var(--ic-neg)"/>
        <circle cx="334" cy="56" r="13" fill="var(--ic-neg)"/>
        <text x="323" y="38" font-size="18" font-weight="600" fill="var(--ic-neg)" text-anchor="middle">−</text>
        <text x="323" y="96" font-size="12" fill="var(--ic-muted)" text-anchor="middle">negative ion</text>
      </svg>
</div>
```

A freed electron quickly loses its extra speed in collisions. In air it then meets an
oxygen molecule that captures it, within about 10 ns and at most a few tens of
micrometres from where it was freed. The charge is not lost. It now sits on a heavy
negative ion, which moves about 200 times slower than the electron did.

So in air, almost all the negative charge crosses the chamber as slow ions. In a gas
with nothing to capture electrons, such as pure nitrogen or argon, the electrons stay
free and reach the anode in well under a millisecond. That is one way the carrier gas
changes the chamber.

**In the model.** Attachment is a reaction that turns electrons into negative ions at
the rate $\nu_a n_e$, with $\nu_a$ the attachment frequency
({class}`~ion_chambers.CarrierReaction`).

## Drift and diffusion

This is the chamber seen from above: the anode rod in the middle at a positive
voltage, the can around it at zero. Ion pairs appear all over the gas. Change the
activity and the voltage, and watch what share of the pairs is collected.

```{raw} html
<div class="ic-stage">
<div class="ic-readouts">
<div class="ic-readout"><span class="ic-label">Collected</span><span class="ic-value" id="ic-sim-collected">0</span></div>
<div class="ic-readout"><span class="ic-label">Recombined</span><span class="ic-value" id="ic-sim-recombined">0</span></div>
<div class="ic-readout"><span class="ic-label">Collection efficiency</span><span class="ic-value" id="ic-sim-efficiency">–</span></div>
</div>
<canvas id="ic-chamber-canvas" aria-label="Animation of positive ions drifting to the can and negative ions drifting to the rod, with some recombining"></canvas>
<div class="ic-controls">
<div class="ic-group" role="group" aria-label="Activity"><span class="ic-label">Activity</span>
<button id="ic-act15" type="button" aria-pressed="false">Low</button>
<button id="ic-act60" type="button" aria-pressed="true">Medium</button>
<button id="ic-act240" type="button" aria-pressed="false">High</button></div>
<div class="ic-group" role="group" aria-label="Voltage"><span class="ic-label">Voltage</span>
<button id="ic-volt200" type="button" aria-pressed="false">200 V</button>
<button id="ic-volt500" type="button" aria-pressed="true">500 V</button>
<button id="ic-volt1500" type="button" aria-pressed="false">1500 V</button></div>
<button id="ic-sim-play" type="button">Pause</button>
</div>
</div>
<div class="ic-legend">
<span><i class="ic-dot" style="background: var(--ic-pos)"></i>positive ion, goes to the can</span>
<span><i class="ic-dot" style="background: var(--ic-neg)"></i>negative ion, goes to the rod</span>
<span><i class="ic-ring"></i>a pair recombines</span>
</div>
```

:::{note}
The animation is slowed down about 100 times: a real ion needs about 44 ms to cross at
500 V. The rod is drawn to scale. The number of ions and the recombination rate are
chosen to make the effect visible, so the percentages show the trend, not the real
values.
:::

The field pushes positive ions outwards to the can (the cathode) and negative ions
inwards to the rod (the anode). An ion is not in free flight. It bumps into a gas
molecule every 70 nm or so, and the field only adds a small steady push between bumps.
The result is a steady drift velocity, proportional to the field:

$$
v_s = -z_s\, \mu_s \nabla \varphi
$$

with $\mu_s$ the mobility, $z_s$ the charge number and $\varphi$ the electric
potential. The field is strongest at the rod, so ions are fast there and slow near the
can.

On top of the drift, the random bumps spread the ions out. This is diffusion, with the
flux $-D_s \nabla n_s$. During one crossing it moves an ion by less than a millimetre,
so it matters little here.

In the model, the diffusion coefficient of a carrier defaults to the Einstein
relation, and the solver raises it where needed to keep the drift stable
({class}`~ion_chambers.DiffusionFloor`, see the
[transport guide](guide/transport.md#stabilisation)). That added diffusion is
numerical, not physical: it makes the error of a current first order in the cell size,
which is why the studies extrapolate from two meshes.

**In the model.** The drift term is FESTIM's `ElectromigrationTerm`, with the velocity
set by the mobility ({class}`~ion_chambers.MobilityElectromigration`). The potential
comes from an {class}`~ion_chambers.ElectrostaticProblem`: see
[the equations together](#the-equations-together).

## Recombination

When a positive and a negative ion come close, the extra electron goes back and both
are neutral molecules again. That ion pair never reaches the electrodes, so its charge
is missing from the reading. In the animation above, a ring marks each one.

- **More activity, more loss.** Twice the activity gives twice as many ions of each
  sign, so four times as many meetings. The share lost grows with the activity. This
  is what bends the response of the chamber at high activity.
- **More voltage, less loss.** Faster ions spend less time in the gas, so fewer are
  in it at once and each has less time to meet a partner.

**In the model.** Recombination is a reaction between two carriers. Positive and
negative ions recombine at the rate $\alpha\, n_P n_N$, positive ions and electrons at
the rate $\beta\, n_P n_e$. In the reference chamber at 1 Ci/m³, 0.7% of the pairs
recombine.

## Collection and the chamber current

Charges that reach an electrode are absorbed by the metal. Their motion through the
gas drives a current in the circuit that connects the rod and the can, and that
current is what the electrometer shows. For the reference chamber at 1 Ci/m³, every
second:

| | |
|---|---:|
| Tritium decays in the gas | 3.69 × 10⁷ |
| × 168 pairs per decay, if nothing were lost | 6.18 × 10⁹ pairs |
| − 1.1% of the energy lost to the walls | 6.12 × 10⁹ pairs |
| − 0.7% of the pairs lost to recombination | 6.08 × 10⁹ pairs |
| × the charge of one electron | **0.97 nA** |

This is the mesh-converged value. A single mesh gives slightly less (0.96 nA in the
walkthrough).

Read from the bottom up, the table turns a current back into decays: divide the
current by the charge of one electron, by the two efficiencies and by 168. The result
is the number of decays per second, and dividing by the gas volume gives the activity
of the gas.

- **It is an average, not a count.** One decay moves about 168 electron charges, far
  too little to see on its own. The chamber measures the steady flow from millions of
  decays per second.
- **It only works if the efficiencies are known.** While they stay constant, the
  current is proportional to the activity and one calibration is enough. Recombination
  grows with the activity, so at some point it is not. Finding that range is what the
  model is for.

**In the model.** An electrode absorbs the carriers the field drives into it (an
outflow condition) and emits none (zero density for the carriers the field pushes
away). Insulators carry no flux. The reading is the current induced on the anode, from
the Shockley–Ramo theorem:

$$
I = e \int \sum_s z_s\, J_s \cdot \left(-\nabla \varphi_w\right) dV,
\qquad \varphi_w = \varphi / V
$$

where $\varphi_w$ is the weighting potential. In steady state this equals the current
through either electrode ({class}`~ion_chambers.ChamberCurrent`).

## How fast it all is

Nothing here moves at the speed of light. The beta is fast, but it is gone almost at
once. Everything after it is slow.

| | Speed | How long it lasts | How far it goes |
|---|---|---|---|
| The beta | 4.4 × 10⁷ m/s (15% of light) | tens of picoseconds | 0.39 mm at the average energy, 4.9 mm at most |
| A free electron in air | 100 to 4 000 m/s with the example mobility | about 10 ns, until it attaches | at most a few tens of micrometres |
| An ion | 0.6 to 19 m/s | about 44 ms to cross the gap | up to 48.5 mm |

The slow step sets the behaviour. Ions stay in the gas for tens of milliseconds, a
billion times longer than the beta that made them. That waiting time is what gives
recombination its chance.

## The equations together

### Charge transport

Each charge carrier $s$ has a number density $n_s$ (m⁻³) that obeys a steady-state
Nernst–Planck balance:

$$
\nabla \cdot J_s = S_s + R_s, \qquad
J_s = -D_s \nabla n_s + n_s v_s, \qquad
v_s = -z_s\, \mu_s \nabla \varphi
$$

with $D_s$ the diffusion coefficient, $S_s$ the ionisation source and $R_s$ the net
production by reactions. This is what {class}`~ion_chambers.ChargeTransportProblem`
solves.

For air, the three balances are

$$
\begin{aligned}
\nabla \cdot J_e &= S - \nu_a n_e - \beta\, n_P n_e \\
\nabla \cdot J_P &= S - \alpha\, n_P n_N - \beta\, n_P n_e \\
\nabla \cdot J_N &= \nu_a n_e - \alpha\, n_P n_N
\end{aligned}
$$

for electrons, positive ions and negative ions.

### Electric potential

The potential is the solution of an {class}`~ion_chambers.ElectrostaticProblem`. It
solves Laplace's equation, $\nabla^2 \varphi = 0$, with $\varphi = V$ on the anode,
$\varphi = 0$ on the cathode and no normal field on the insulators. The charge in the
gas is not fed back into the potential.

To check that assumption, the potential $\varphi_1$ created by the converged charge
density is computed afterwards from

$$
-\nabla \cdot (\varepsilon_0 \nabla \varphi_1) = e \sum_s z_s n_s,
\qquad \varphi_1 = 0 \text{ on the electrodes}
$$

and the ratio $\lVert \nabla\varphi_1 \rVert / \lVert \nabla\varphi \rVert$ is
reported by the export {class}`~ion_chambers.SpaceChargeField`.

## What the model leaves out

- **Space-charge feedback.** The field is fixed. The distortion is reported but not
  applied.
- **Gas flow.** The gas is still. In a flow-through chamber the flow is two to four
  orders of magnitude slower than the ion drift.
- **Columnar recombination** along individual beta tracks. Only volume recombination
  between uniformly mixed carriers is modelled.
- **Tritium on the walls**, which causes memory effects and background drift. The
  [walkthrough](walkthrough.md) sketches how FESTIM could be coupled for this.
- **Time dependence.** Only the steady state is solved.

```{raw} html
<script src="_static/mechanisms.js"></script>
```
