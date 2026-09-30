# neutrino_oscillations

A small, self-contained Python package (`nuosc`) to compute three-flavour
neutrino oscillation probabilities **in vacuum** and **in constant-density
matter**, for all nine channels $P(\nu_\alpha \to \nu_\beta)$, plus a notebook
with example plots. Written to be read and modified by bachelor students.

## Contents

```
neutrino_oscillations/
├── nuosc/
│   ├── constants.py     # physical constants, unit conventions
│   ├── parameters.py    # NuFIT 6.0 best-fit values for NO and IO
│   ├── experiments.py   # experiment configs: T2K, NOvA, DUNE, JUNO (+ flux, mass, POT)
│   ├── formulas.py      # vacuum formula, NuFast (matter), exact numerical check
│   ├── oscillator.py    # NeutrinoOscillator class (the thing you normally use)
│   ├── fluxes.py        # loaders for the far-detector fluxes
│   ├── rates.py         # cross section, event rates, toy energy smearing
│   ├── data.py          # NOvA 2024 and T2K Run 1-9 far-detector data
│   ├── dune_globes.py   # DUNE TDR simulation (official GLoBES configuration in Python)
│   ├── spectra.py       # load / oscillate the far-detector spectra below
│   ├── plotting.py      # optional set_style(True/False): latex_style -> SciencePlots -> default
│   └── paths.py         # finds the data folders
├── fluxes/              # public unoscillated fluxes (see fluxes/README.md)
│   ├── DUNE/  ├── NOvA/  └── T2K/
├── far_detector_spectra/  # T2K, NOvA, DUNE, HyperK: FD spectra ready to multiply by P (one CSV each)
├── tools/build_fd_spectra.py  # rebuilds far_detector_spectra/
├── data/                # public data and configurations (see data/README.md)
│   ├── NOvA_2024_FD_release/  ├── T2K_OA2019_release/  ├── DUNE_GLoBES_2103.04797/  └── HyperK_DR_1805.04163/
├── pyproject.toml               # pip install -e .
├── GUIDE.md                      # detailed user guide: nu/nubar, NO/IO, matter, parameters, experiments
├── tutorial_nuosc.ipynb          # START HERE: how to use the code (probabilities only) + exercises
├── tutorial_exercises_results.ipynb  # results of the two tutorial exercises
├── solutions_matter_biprobability.ipynb  # beyond the tutorial: matter effects, bi-probability, bi-event plots, delta_CP
├── oscillation_plots.ipynb       # probabilities
├── realistic_event_rates.ipynb   # fluxes x probabilities x cross section = events; NOvA & T2K data
├── dune_tdr_simulation.ipynb     # DUNE with the official TDR GLoBES configuration
├── delta_cp_LE_interference.ipynb  # why P(numubar->nuebar) vs L/E looks "strange": delta_CP
├── fd_spectra_quickstart.ipynb   # use far_detector_spectra/ for all four experiments
└── requirements.txt              # same dependencies as pyproject.toml
```

## Installation

From this folder:

```bash
pip install -e .
```

This installs the `nuosc` package and everything the notebooks need (numpy, matplotlib,
SciencePlots, jupyter). `-e` (editable) keeps the package linked to this folder, so the data
folders (`fluxes/`, `data/`, `far_detector_spectra/`) are found automatically; with a plain
`pip install .` run the notebooks from this folder or set `NUOSC_ROOT` to its path.

**Plot style is optional.** Every notebook starts with `SET_STYLE = True` and `c = set_style(SET_STYLE)`:

* `True`: the project style — the `latex_style` matplotlib style if installed, otherwise SciencePlots,
  otherwise the matplotlib default, with the IBM colour-blind-safe colour cycle;
* `False`: nothing is changed, your own matplotlib settings are used.

`c` is always the list of colours of the active colour cycle. Another style can be tried first with
`set_style(True, style="my_style")`, and `colors=None` keeps that style's own colours.

**Start here:** `tutorial_nuosc.ipynb` (how to use the code, with exercises), results in
`tutorial_exercises_results.ipynb`. The detailed reference is [`GUIDE.md`](GUIDE.md) (how to switch
$\nu/\bar\nu$, NO/IO, vacuum/matter, change parameters and experiments). Bi-probability plots and the
$\delta_{CP}$ measurement are in `solutions_matter_biprobability.ipynb`.

## Units

| quantity        | unit    |
|-----------------|---------|
| energy `E`      | GeV (reactor: `E = E_MeV * 1e-3`) |
| baseline `L`    | km      |
| `dm2_21`, `dm2_31` | eV²  |
| density `rho`   | g/cm³   |
| `delta_cp`      | degrees |

## Quick start

```python
import numpy as np
from nuosc import NeutrinoOscillator

osc = NeutrinoOscillator(ordering="NO", baseline=1300, rho=2.848)  # DUNE-like
print(osc)                                   # current parameters

E = np.linspace(0.5, 5, 500)                 # GeV
P_vac = osc.vacuum("mu", "e", E)             # P(nu_mu -> nu_e), vacuum
P_mat = osc.matter("mu", "e", E)             # P(nu_mu -> nu_e), matter (NuFast)
P_all = osc.matter(E=E)                      # all channels, shape (3, 3, 500)

osc.ordering = "IO"                          # switch ordering (reloads NuFIT IO)
osc.set_params(delta_cp=270, s2_23=0.5)      # change any parameter
osc.baseline, osc.rho = 295, 2.6             # change experiment
osc.antineutrino = True                      # antineutrinos
osc.reset()                                  # back to NuFIT defaults
```

### Experiments

Experiment setups (baseline, density, energy range, beam mode, main channel)
are stored in `nuosc/experiments.py`:

```python
osc = NeutrinoOscillator.from_experiment("NOvA", ordering="IO", antineutrino=True)
osc.matter("mu", "e")                        # on the experiment's energy grid
```

Available: `T2K`, `NOvA`, `DUNE`, `HyperK`, `JUNO`. Add your own by adding an entry to
`EXPERIMENTS`. In the notebook, change `EXPERIMENT` / `COMPARE` in the
configuration cell and re-run.

`E` and `L` broadcast against each other, so you can pass 2D meshgrids to
make oscillograms.

### Fluxes and event rates

```python
from nuosc import load_flux, expected_events, group_rates, integrate

flux = load_flux("DUNE", "FHC")              # nu / cm^2 / GeV / POT at the FD
E, rates = expected_events("DUNE", "FHC", osc=NeutrinoOscillator("NO"))
samples = group_rates(rates, "FHC")          # mu_like, e_like, e_signal, e_beam, ...
integrate(E, samples["e_like"])              # total nu_e-like CC events
```

Event rates use a linear CC cross section and efficiency 1; see `nuosc/rates.py`.
The NOvA far-detector flux is approximated from the public near-detector flux
(see `fluxes/README.md`); for NOvA the notebook also compares with the
real far-detector data of the 2024 release (`nuosc/data.py`).

### DUNE TDR simulation

```python
from nuosc.dune_globes import DUNEGlobes
dune = DUNEGlobes()                          # reads data/DUNE_GLoBES_2103.04797/dune_globes
s = dune.sample(NeutrinoOscillator("NO"), "nue_app")
s["signal"].sum(), s["background"].sum()     # selected events in reconstructed-energy bins
```

Flux, GENIE cross sections on argon, smearing matrices and efficiencies are those of the DUNE
GLoBES configuration (arXiv:2103.04797); the notebook reproduces the spectra of that paper.

### Far-detector spectra (the easy way)

```python
from nuosc import predict, NeutrinoOscillator
s = predict("HyperK", NeutrinoOscillator("NO"))   # also "T2K", "NOvA", "DUNE"
s["FHC"]["e_like"].sum(), s["RHC"]["mu_like"].sum()
```

Each `far_detector_spectra/<EXP>_FD_spectra.csv` lists selected events per energy bin for P = 1,
split by flavour; multiply by the probabilities and add (recipe in the file header and in
`far_detector_spectra/README.md`). The files can be used without `nuosc`.

## Formulas

* **Vacuum** – the standard PMNS expression
  $P_{\alpha\beta} = \delta_{\alpha\beta} - 4\sum_{i>j}\mathrm{Re}(U^*_{\alpha i}U_{\beta i}U_{\alpha j}U^*_{\beta j})\sin^2\Delta_{ij} + 2\sum_{i>j}\mathrm{Im}(\dots)\sin 2\Delta_{ij}$,
  with $\Delta_{ij} = 1.267\,\Delta m^2_{ij}[\mathrm{eV^2}]\,L[\mathrm{km}]/E[\mathrm{GeV}]$.
* **Matter** – the NuFast algorithm, P. B. Denton and S. J. Parke,
  [arXiv:2405.02400](https://arxiv.org/abs/2405.02400). `N_Newton = 0` is
  already accurate to ~1e-5 or better; `N_Newton = 1` gives ~1e-9.
* **Exact** – `method="exact"` diagonalises the Hamiltonian numerically;
  slower, useful as a cross-check.

Antineutrinos: $U \to U^*$ ($\delta_{CP}\to-\delta_{CP}$) and $A \to -A$.

Default parameters: NuFIT 6.0 (IC24 with SK atmospheric data),
[arXiv:2410.05380](https://arxiv.org/abs/2410.05380). Edit `nuosc/parameters.py`
to use other values.
