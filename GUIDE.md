# `nuosc` user guide

This guide explains how to use the code in detail: what every option does, how to switch between neutrinos and
antineutrinos, between normal and inverted ordering, between vacuum and matter, and how to change parameters
and experiments. Each section has the physics behind it in a few lines, the code, and common pitfalls.

For a hands-on introduction, run `tutorial_nuosc.ipynb`; this file is the reference to come back to.

**Contents**

1. [Installation and first lines](#1-installation-and-first-lines)
2. [Units](#2-units)
3. [The `NeutrinoOscillator` object](#3-the-neutrinooscillator-object)
4. [Channels and the probability matrix](#4-channels-and-the-probability-matrix)
5. [Vacuum or matter](#5-vacuum-or-matter)
6. [Neutrinos or antineutrinos](#6-neutrinos-or-antineutrinos)
7. [Normal or inverted ordering](#7-normal-or-inverted-ordering)
8. [Changing the oscillation parameters](#8-changing-the-oscillation-parameters)
9. [Baseline, energy and density](#9-baseline-energy-and-density)
10. [Experiment configurations](#10-experiment-configurations)
11. [Scans: energy, $L/E$, parameters](#11-scans-energy-le-parameters)
12. [Keeping several configurations side by side](#12-keeping-several-configurations-side-by-side)
13. [Plot style](#13-plot-style)
14. [Checks you can do yourself](#14-checks-you-can-do-yourself)
15. [Common mistakes and FAQ](#15-common-mistakes-and-faq)
16. [Where the code lives](#16-where-the-code-lives)

---

## 1. Installation and first lines

From the project folder:

```bash
pip install -e .
```

Then, in Python or in a notebook:

```python
import numpy as np
from nuosc import NeutrinoOscillator

osc = NeutrinoOscillator()          # NO, L = 1300 km, rho = 2.8 g/cm^3, neutrinos, NuFIT 6.0
print(osc)                          # shows all parameters
osc.matter("mu", "e", E=2.5)        # P(nu_mu -> nu_e) at 2.5 GeV, in matter
```

## 2. Units

| quantity | unit | example |
|---|---|---|
| energy `E` | GeV | `E=0.6` (not 600 MeV!) |
| baseline `L`, `baseline` | km | `295`, `810`, `1300` |
| `dm2_21`, `dm2_31`, `dm2_32` | eV² | `2.5e-3` |
| density `rho` | g/cm³ | `2.8` (Earth's crust), `0` = vacuum |
| `delta_cp` | degrees | `270` (not $3\pi/2$!) |
| mixing angles | $\sin^2\theta_{ij}$ | `s2_23=0.5` |

## 3. The `NeutrinoOscillator` object

Everything goes through one object that stores the **settings** (ordering, baseline, density, beam type)
and the **six oscillation parameters**.

```python
osc = NeutrinoOscillator(
    ordering="NO",         # "NO" or "IO": loads the NuFIT 6.0 best fit for that ordering
    baseline=1300,         # L [km]
    energy=None,           # default energy grid [GeV]; None -> np.linspace(0.2, 10, 2000)
    rho=2.8,               # matter density [g/cm^3]
    Ye=0.5,                # electrons per nucleon (0.5 in the Earth)
    antineutrino=False,    # False: neutrinos, True: antineutrinos
    N_Newton=0,            # NuFast accuracy (0 is ~1e-5 already, 1 is ~1e-9)
    delta_cp=270,          # optional: any parameter can be overridden here
)
```

Everything can be read and changed later as an attribute:

```python
osc.baseline, osc.rho, osc.antineutrino, osc.ordering, osc.energy
osc.s2_12, osc.s2_13, osc.s2_23, osc.delta_cp, osc.dm2_21, osc.dm2_31
osc.dm2_32                  # derived: dm2_31 - dm2_21
osc.theta_12, osc.theta_13, osc.theta_23   # angles in degrees (read only)
osc.params                  # dict with the six parameters
osc.U                       # 3x3 complex PMNS matrix
```

## 4. Channels and the probability matrix

Flavours are called `"e"`, `"mu"`, `"tau"` (or `0`, `1`, `2`). The **first argument is the initial flavour**, the
second the final one.

```python
osc.matter("mu", "e", E)     # P(nu_mu -> nu_e)   appearance
osc.matter("mu", "mu", E)    # P(nu_mu -> nu_mu)  disappearance (survival)
osc.matter("e", "e", E)      # P(nu_e -> nu_e)

P = osc.matter(E=E)          # all nine at once: P[initial, final, ...]
P[1, 0]                      # = P(nu_mu -> nu_e)
P.sum(axis=1)                # = 1 for every initial flavour (unitarity)
```

`E` can be a number or an array; the result has the same shape.

## 5. Vacuum or matter

**Physics.** On the way through the Earth, $\nu_e$ scatter coherently on electrons (charged current).
This adds a potential $A = 2\sqrt2\,G_F N_e E$ to the $\nu_e$ entry of the Hamiltonian. It changes the effective
mixing angles and mass splittings. The effect grows with $E$ (and so, at the oscillation maximum, with $L$):
it is small at T2K (295 km), sizeable at NOvA (810 km), large at DUNE (1300 km).

**Code.** Three equivalent ways:

```python
osc.vacuum("mu", "e", E)                       # vacuum formula
osc.matter("mu", "e", E)                       # constant density osc.rho (NuFast)
osc.probability("mu", "e", E, matter=False)    # same as vacuum
osc.probability("mu", "e", E, matter=True)     # same as matter

osc.matter("mu", "e", E, method="exact")       # numerical diagonalisation (slow, cross-check)
```

Setting `osc.rho = 0` also gives the vacuum result through the matter formula (a good check).

The size of the effect is measured by $A/\Delta m^2_{31}$:

```python
from nuosc import matter_potential
A = matter_potential(E=2.5, rho=2.848)         # eV^2, sign flips for antineutrinos
print(A / osc.dm2_31)                          # ~0.2 at DUNE
```

## 6. Neutrinos or antineutrinos

**Physics.** Going from $\nu$ to $\bar\nu$ changes two things:

1. $U \to U^*$, i.e. $\delta_{CP} \to -\delta_{CP}$. This is **genuine CP violation**: it exists also in vacuum
   and it is what we want to measure.
2. $A \to -A$: the matter potential changes sign. This is a "fake" CP asymmetry, due to the Earth being made of
   matter and not antimatter. It exists even if $\delta_{CP}=0$.

So in vacuum with $\delta_{CP}=0$ or $180°$, $P(\nu_\mu\to\nu_e) = P(\bar\nu_\mu\to\bar\nu_e)$ exactly; in matter
they differ anyway. Separating the two effects is the reason experiments need a long baseline **and** both beam
modes.

**Code.** The switch is one attribute. You do **not** change `delta_cp` or `rho` yourself: the code applies
both changes internally.

```python
osc.antineutrino = True       # from now on: P(nubar_alpha -> nubar_beta)
osc.matter("mu", "e", E)      # = P(nubar_mu -> nubar_e)
osc.antineutrino = False      # back to neutrinos
```

or directly at creation:

```python
nu    = NeutrinoOscillator("NO", baseline=810, rho=2.84)
nubar = NeutrinoOscillator("NO", baseline=810, rho=2.84, antineutrino=True)
asym  = (nu.matter("mu", "e", E) - nubar.matter("mu", "e", E)) / (nu.matter("mu", "e", E) + nubar.matter("mu", "e", E))
```

`print(osc)` tells you which one you are using (`neutrino` / `antineutrino`).

In the experiment configurations `antineutrino` is the default beam mode (`False`, i.e. the neutrino beam, FHC).
Override it with
`NeutrinoOscillator.from_experiment("NOvA", antineutrino=True)` (RHC beam).

## 7. Normal or inverted ordering

**Physics.** We know $|\Delta m^2_{31}|$ but not its sign:

* **NO** (normal ordering): $m_1 < m_2 < m_3$, $\Delta m^2_{31} > 0$;
* **IO** (inverted ordering): $m_3 < m_1 < m_2$, $\Delta m^2_{31} < 0$.

In vacuum the sign enters only through the phases and is hard to see. In matter it matters a lot: the matter
effect **enhances** $\nu_\mu\to\nu_e$ for NO and **suppresses** it for IO, and the opposite for antineutrinos.

**Code.**

```python
osc = NeutrinoOscillator(ordering="IO")   # at creation
osc.ordering = "IO"                       # or later
osc.ordering                              # -> "IO"
```

**Important:** setting the ordering **reloads all six parameters** with the NuFIT 6.0 best fit for that ordering.
NuFIT's IO best fit is not just "NO with the sign flipped":

| | NO | IO |
|---|---|---|
| $\sin^2\theta_{12}$ | 0.308 | 0.308 |
| $\sin^2\theta_{13}$ | 0.02215 | 0.02231 |
| $\sin^2\theta_{23}$ | 0.470 | 0.550 |
| $\delta_{CP}$ | 212° | 274° |
| $\Delta m^2_{21}$ | 7.49e-5 eV² | 7.49e-5 eV² |
| $\Delta m^2_{31}$ | +2.513e-3 eV² | −2.409e-3 eV² ($\Delta m^2_{32} = -2.484$e-3) |

So:

* change the ordering **first**, then your own parameters:
  ```python
  osc.ordering = "IO"
  osc.set_params(delta_cp=270)      # correct order
  ```
* if you want to isolate the effect of the ordering alone (same angles and phase, only the sign of the
  splitting), copy the parameters and flip the sign yourself:
  ```python
  no = NeutrinoOscillator("NO")
  io = NeutrinoOscillator("IO")
  io.set_params(**{k: v for k, v in no.params.items() if k != "dm2_31"})   # same angles, delta, dm2_21
  io.set_params(dm2_32=-no.dm2_32)                                          # same |dm2_32|, opposite sign
  ```
* `set_params` prints a warning if the sign of `dm2_31` does not match the ordering (e.g. you set a positive
  `dm2_31` while in IO). The probability is still computed with the numbers you gave.

The NuFIT values live in `nuosc/parameters.py`; edit them there if you want other defaults.

## 8. Changing the oscillation parameters

```python
osc.set_params(delta_cp=270)                   # one parameter
osc.set_params(s2_23=0.5, dm2_31=2.5e-3)       # several
osc.set_params(dm2_32=2.43e-3)                 # dm2_32 is accepted too (converted to dm2_31)
osc.delta_cp = 0                               # direct assignment also works
osc.reset()                                    # back to NuFIT for the current ordering
```

Valid names: `s2_12`, `s2_13`, `s2_23`, `delta_cp`, `dm2_21`, `dm2_31`, `dm2_32`. A typo raises a `KeyError`
listing the valid names.

What each parameter mainly does in the accelerator channels:

| parameter | effect |
|---|---|
| `s2_13` | overall size of $P(\nu_\mu\to\nu_e)$ ($\propto\sin^2 2\theta_{13}$) |
| `s2_23` | depth of the $P(\nu_\mu\to\nu_\mu)$ dip ($\sin^2 2\theta_{23}$) and size of $P(\nu_\mu\to\nu_e)$ ($\sin^2\theta_{23}$, octant) |
| `dm2_31` / `dm2_32` | position in energy of the maxima and minima; its sign is the ordering |
| `delta_cp` | $\nu$ vs $\bar\nu$ difference in $P(\nu_\mu\to\nu_e)$, and a shift of the peak; almost no effect on disappearance |
| `s2_12`, `dm2_21` | solar parameters: small effects at these baselines (they enter the $\delta_{CP}$ interference term) |

## 9. Baseline, energy and density

```python
osc.baseline = 810                    # km
osc.rho = 2.84                        # g/cm^3
osc.energy = np.linspace(0.5, 5, 500) # default grid used when E is not given

osc.matter("mu", "e")                 # uses osc.energy and osc.baseline
osc.matter("mu", "e", E=2.0)          # one energy
osc.matter("mu", "e", E=2.0, L=500)   # another baseline, only for this call
```

`E` and `L` broadcast like numpy arrays, so `E` and `L` can be arrays of the same shape, or a 2D grid from
`np.meshgrid`.

## 10. Experiment configurations

`nuosc/experiments.py` stores baseline, density, typical energy and energy range of the long-baseline experiments:

| name | $L$ [km] | $\rho$ [g/cm³] | `E_peak` [GeV] | `E_range` [GeV] | default beam |
|---|---|---|---|---|---|
| `T2K` | 295 | 2.6 | 0.6 | 0.1–3 | $\nu$ |
| `HyperK` | 295 | 2.6 | 0.6 | 0.1–3 | $\nu$ |
| `NOvA` | 810 | 2.84 | 2.0 | 0.3–6 | $\nu$ |
| `DUNE` | 1300 | 2.848 | 2.5 | 0.3–8 | $\nu$ |

```python
from nuosc import get_experiment

osc = NeutrinoOscillator.from_experiment("NOvA")                       # L, rho, energy grid, beam set
osc = NeutrinoOscillator.from_experiment("NOvA", ordering="IO", antineutrino=True, delta_cp=90)

cfg = get_experiment("DUNE")          # the dictionary itself
cfg["baseline"], cfg["rho"], cfg["E_peak"], cfg["E_range"]
osc.matter("mu", "e", E=cfg["E_peak"])  # probability at the typical energy

for name in ["T2K", "HyperK", "NOvA", "DUNE"]:   # loop over experiments
    osc = NeutrinoOscillator.from_experiment(name)
```

`E_peak` is the typical energy of the beam, close to the mean energy of the appeared $\nu_e$ events.
To add an experiment, add an entry to `EXPERIMENTS` with at least `baseline`, `rho`, `E_peak`, `E_range`,
`antineutrino`, `channel`.

## 11. Scans: energy, $L/E$, parameters

**Energy at fixed baseline**

```python
E = np.linspace(0.3, 8, 1000)
P = osc.matter("mu", "e", E)                  # array of 1000 values
```

**$L/E$.** In vacuum $P$ depends only on $L/E$. Two ways:

```python
LE = osc.baseline / E                         # fixed L, vary E
plt.plot(LE, osc.vacuum("mu", "mu", E))

LE = np.linspace(10, 2000, 1000)              # or fix E and vary L = LE * E
plt.plot(LE, osc.matter("mu", "mu", E=2.0, L=LE * 2.0))
```

In matter the two are different (the matter term depends on $E$ alone), so say which one you use.

**A parameter**

```python
deltas = np.arange(0, 361, 5)
P = []
for d in deltas:
    osc.set_params(delta_cp=d)
    P.append(osc.matter("mu", "e", E=2.5))
osc.reset()                                   # do not forget!
```

**Position of the first maximum.** Useful in exercises:

```python
E = np.linspace(1, 5, 4000)
P = osc.matter("mu", "e", E)
E_max = E[np.argmax(P)]
```

Choose the range so that it contains only the maximum you want.

## 12. Keeping several configurations side by side

Changing attributes in place is easy but can leave the object in a state you forgot about. For comparisons it is
often clearer to create one object per configuration:

```python
configs = {
    ("NO", False): NeutrinoOscillator.from_experiment("DUNE", ordering="NO"),
    ("NO", True):  NeutrinoOscillator.from_experiment("DUNE", ordering="NO", antineutrino=True),
    ("IO", False): NeutrinoOscillator.from_experiment("DUNE", ordering="IO"),
    ("IO", True):  NeutrinoOscillator.from_experiment("DUNE", ordering="IO", antineutrino=True),
}
for (o, anti), osc in configs.items():
    plt.plot(osc.energy, osc.matter("mu", "e"), ls="--" if anti else "-", label=f"{o}, {'nubar' if anti else 'nu'}")
```

## 13. Plot style

```python
from nuosc.plotting import set_style
SET_STYLE = True
c = set_style(SET_STYLE)
```

* `True`: `latex_style` if you have it, otherwise SciencePlots, otherwise matplotlib's default, with the IBM
  colour-blind-safe colour cycle.
* `False`: nothing is changed; your own matplotlib settings are used.
* `c` is the list of colours (`c[0]`, `c[1]`, ...) in both cases.
* `set_style(True, style="my_style")` tries your style first; `colors=None` keeps that style's own colours.

## 14. Checks you can do yourself

```python
P = osc.matter(E=E);  np.allclose(P.sum(axis=1), 1)               # unitarity
osc.rho = 0;  np.allclose(osc.matter("mu", "e", E), osc.vacuum("mu", "e", E))
np.abs(osc.matter("mu", "e", E) - osc.matter("mu", "e", E, method="exact")).max()   # ~1e-5
# vacuum, delta_CP = 0: nu and nubar identical
```

## 15. Common mistakes and FAQ

| mistake | what happens | fix |
|---|---|---|
| energy in MeV | `E=600` is 600 GeV | GeV: `E=0.6` |
| `delta_cp` in radians | `delta_cp=np.pi` is 3.14° | degrees |
| ordering changed after the parameters | your parameters are overwritten by NuFIT | ordering first |
| `antineutrino=True` left on | later cells are for $\bar\nu$ | set it back, or use two objects |
| parameters changed in a loop | following cells use the last value | `osc.reset()` |
| `P[1, 0]` vs `P[0, 1]` | $\nu_\mu\to\nu_e$ vs $\nu_e\to\nu_\mu$ | first index = initial flavour |
| flipping `delta_cp` by hand for $\bar\nu$ | the sign is flipped twice | only set `antineutrino=True` |

**Why are T2K and Hyper-K identical?** Same baseline, same beam, same density: the probabilities are the same;
they differ in detector size (number of events), not in oscillation physics.

**Why does $P(\nu_\mu\to\nu_\mu)$ barely depend on $\delta_{CP}$ and on matter?** Its leading term is
$\sin^2 2\theta_{23}\sin^2\Delta_{31}$; $\delta_{CP}$ enters only through small $\cos\delta_{CP}$ terms, and the matter
potential acts on $\nu_e$, which is only a small part of this channel.

**Why does $P(\nu_e\to\nu_e)$ not depend on $\delta_{CP}$ at all?** $\delta_{CP}$ can be moved into the
$\mu$–$\tau$ sector by a rephasing that does not touch $\nu_e$, both in vacuum and in constant-density matter.

**Is NuFast exact?** For constant density it is accurate to ~$10^{-5}$ with `N_Newton=0` and ~$10^{-9}$ with
`N_Newton=1`; `method="exact"` diagonalises the Hamiltonian numerically.

**Real Earth density?** The code uses one constant density along the path, a very good approximation for these
baselines (the beam only crosses the crust).

## 16. Where the code lives

| file | content |
|---|---|
| `nuosc/oscillator.py` | `NeutrinoOscillator` (what you normally use) |
| `nuosc/formulas.py` | PMNS matrix, vacuum formula, NuFast, exact diagonalisation |
| `nuosc/parameters.py` | NuFIT 6.0 default values for NO and IO |
| `nuosc/experiments.py` | experiment configurations |
| `nuosc/constants.py` | physical constants and unit conversions |
| `nuosc/plotting.py` | `set_style` |

Notebooks: `tutorial_nuosc.ipynb` (start here), `tutorial_exercises_results.ipynb` (results of the tutorial
exercises), then the more advanced ones listed in `README.md`.
