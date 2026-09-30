"""
Experiment configurations
=========================

One place to store the setup of each experiment. Every entry has

    baseline     : L [km]
    rho          : average matter density along the path [g/cm^3]
    E_peak       : typical (flux-peak) energy [GeV]
    E_range      : (E_min, E_max) energy range for plots [GeV]
    antineutrino : default beam mode (False = neutrinos)
    channel      : main channel (alpha, beta), e.g. ("mu", "e") = nu_mu -> nu_e

Accelerator experiments also have (used by nuosc.fluxes and nuosc.rates)

    flux         : format + file names (relative to the fluxes/ folder) of the
                   unoscillated far-detector flux for FHC (neutrino) and
                   RHC (antineutrino) horn current
    mass_kt      : far-detector target mass [kt]
    pot          : exposure in protons on target, for FHC and RHC
    resolution   : toy fractional energy resolution sigma_E / E

The exposures and masses are round numbers close to the published ones
(T2K Run 1-10, NOvA 2024 analysis, DUNE TDR ~3.5 yr per mode at 1.2 MW,
Hyper-K Design Report 10 yr at 1.3 MW);
the resolutions are simple toy values. Change them freely.

To add an experiment just add a new entry to EXPERIMENTS, then use
    NeutrinoOscillator.from_experiment("MyExperiment")
"""

EXPERIMENTS = {
    "T2K": dict(
        baseline=295.0, rho=2.6,
        E_peak=0.6, E_range=(0.1, 3.0),
        antineutrino=False, channel=("mu", "e"),
        flux=dict(format="t2k_sk",
                  FHC="T2K/t2kflux_2020_public_release/t2kflux_2020_plus250kA_nominal_sk.txt",
                  RHC="T2K/t2kflux_2020_public_release/t2kflux_2020_minus250kA_nominal_sk.txt"),
        mass_kt=22.5,                          # Super-K fiducial volume (water)
        pot=dict(FHC=2.0e21, RHC=1.6e21),
        resolution=0.10,
    ),
    "NOvA": dict(
        baseline=810.0, rho=2.84,
        E_peak=2.0, E_range=(0.3, 6.0),
        antineutrino=False, channel=("mu", "e"),
        flux=dict(format="nova_nd", FHC="NOvA", RHC="NOvA",
                  kwargs=dict(L_ND=1.0, L_FD=810.0)),
        mass_kt=14.0,                          # FD total mass (liquid scintillator)
        pot=dict(FHC=26.61e20, RHC=12.5e20),    # 2024 analysis
        resolution=0.09,
    ),
    "DUNE": dict(
        baseline=1300.0, rho=2.848,
        E_peak=2.5, E_range=(0.3, 8.0),
        antineutrino=False, channel=("mu", "e"),
        flux=dict(format="dune_globes",
                  FHC="DUNE/DUNE_FD_FHC_flux_globes.txt",
                  RHC="DUNE/DUNE_FD_RHC_flux_globes.txt"),
        mass_kt=40.0,                          # 4 x 10 kt LAr far-detector modules
        pot=dict(FHC=3.85e21, RHC=3.85e21),    # 3.5 yr x 1.1e21 POT/yr per mode
        resolution=0.15,
    ),
    "HyperK": dict(                     # same J-PARC beam and off-axis angle as T2K
        baseline=295.0, rho=2.6,
        E_peak=0.6, E_range=(0.1, 3.0),
        antineutrino=False, channel=("mu", "e"),
        flux=dict(format="t2k_sk",
                  FHC="T2K/t2kflux_2020_public_release/t2kflux_2020_plus250kA_nominal_sk.txt",
                  RHC="T2K/t2kflux_2020_public_release/t2kflux_2020_minus250kA_nominal_sk.txt"),
        mass_kt=187.0,                         # one tank, fiducial (Design Report, arXiv:1805.04163)
        pot=dict(FHC=0.675e22, RHC=2.025e22),  # 2.7e22 POT (1.3 MW x 10 yr), nu:nubar = 1:3
        resolution=0.10,
    ),
    "JUNO": dict(                       # reactor: E in GeV (1.8 MeV = 1.8e-3 GeV)
        baseline=52.5, rho=2.45,
        E_peak=4.0e-3, E_range=(1.8e-3, 8.0e-3),
        antineutrino=True, channel=("e", "e"),
    ),
}

# Experiments with a far-detector flux available
ACCELERATORS = [name for name, cfg in EXPERIMENTS.items() if "flux" in cfg]


def get_experiment(name):
    """Return a copy of the configuration of experiment `name` (case-insensitive)."""
    for key, cfg in EXPERIMENTS.items():
        if key.lower() == name.lower():
            return dict(cfg, name=key)
    raise KeyError(f"Unknown experiment {name!r}. Available: {list(EXPERIMENTS)}")
