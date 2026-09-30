"""
Event rates at the far detector
===============================

The expected number of charged-current (CC) events of flavour beta, produced
by the flux of flavour alpha, per unit (true) neutrino energy is

    dN/dE = Phi_alpha(E) * P(alpha -> beta; E) * sigma_CC,beta(E) * N_target * POT * eff

  Phi_alpha : unoscillated flux at the far detector  [nu / cm^2 / GeV / POT]
  P         : oscillation probability (vacuum or matter, from nuosc)
  sigma_CC  : CC cross section per nucleon           [cm^2]
  N_target  : number of nucleons in the fiducial mass
  POT       : exposure, protons on target
  eff       : detection efficiency (default 1 -> "true" interaction rate)

Cross section
-------------
We use the simplest approximation, valid in the deep-inelastic regime and
good to ~20-30% down to ~0.5 GeV for an isoscalar target:

    sigma_CC(nu)    ~ 0.70e-38 cm^2 * E[GeV]
    sigma_CC(nubar) ~ 0.35e-38 cm^2 * E[GeV]

(see the PDG review "Neutrino cross section measurements"). The same value
is used for nu_e and nu_mu; nu_tau CC events are not considered.
"""

import numpy as np

from .constants import N_A
from .experiments import get_experiment
from .fluxes import FLAVOUR_INFO, load_flux
from .formulas import prob_vacuum, prob_matter_nufast

SIGMA_OVER_E_NU = 0.70e-38      # cm^2 / GeV / nucleon
SIGMA_OVER_E_NUBAR = 0.35e-38   # cm^2 / GeV / nucleon

NUCLEONS_PER_KT = 1e9 * N_A     # 1 kt = 1e9 g, ~N_A nucleons per gram


def cross_section_cc(E, antineutrino=False):
    """Approximate CC inclusive cross section per nucleon [cm^2], E in GeV."""
    E = np.asarray(E, float)
    return (SIGMA_OVER_E_NUBAR if antineutrino else SIGMA_OVER_E_NU) * E


def n_nucleons(mass_kt):
    """Number of target nucleons in `mass_kt` kilotons."""
    return mass_kt * NUCLEONS_PER_KT


def final_flavour(initial, final):
    """'numu' + 'e' -> 'nue'; keeps the neutrino / antineutrino nature."""
    bar = "bar" if FLAVOUR_INFO[initial][1] else ""
    return {"e": "nue", "mu": "numu", "tau": "nutau"}[final] + bar


def probability(osc, initial, final, E, L, rho, matter=True):
    """
    P(initial -> final) for a flux flavour name (e.g. 'numubar') and a final
    flavour ('e', 'mu', 'tau'), using the oscillation parameters of `osc`.
    `osc = None` means no oscillations (P = delta).
    """
    a, anti = FLAVOUR_INFO[initial]
    b = {"e": 0, "mu": 1, "tau": 2}[final]
    if osc is None:
        return np.full_like(np.asarray(E, float), 1.0 if a == b else 0.0)
    if matter:
        P = prob_matter_nufast(E, L, **osc.params, rho=rho, Ye=osc.Ye,
                               antineutrino=anti, N_Newton=osc.N_Newton)
    else:
        P = prob_vacuum(E, L, **osc.params, antineutrino=anti)
    return P[a, b]


def expected_events(experiment, mode="FHC", osc=None, E=None, matter=True,
                    efficiency=1.0, pot=None, mass_kt=None, flux=None):
    """
    Differential event rates dN/dE [events / GeV] at the far detector.

    Parameters
    ----------
    experiment : name in nuosc.experiments (e.g. 'DUNE')
    mode       : 'FHC' (neutrino beam) or 'RHC' (antineutrino beam)
    osc        : NeutrinoOscillator (its parameters/ordering are used; baseline
                 and density are taken from the experiment). None = no oscillations
    E          : true-energy grid [GeV]; default: fine grid over E_range
    matter     : include matter effects
    efficiency : scalar or dict {channel: eff}
    pot, mass_kt, flux : override the experiment defaults

    Returns
    -------
    E, rates : rates is a dict {'numu->nue': dN/dE, ...} for all combinations
               of the 4 flux flavours (numu, numubar, nue, nuebar) and the
               2 detected flavours (e, mu).
    """
    cfg = get_experiment(experiment)
    flux = load_flux(experiment, mode) if flux is None else flux
    pot = cfg["pot"][mode] if pot is None else pot
    mass_kt = cfg["mass_kt"] if mass_kt is None else mass_kt
    if E is None:
        E = np.linspace(cfg["E_range"][0], cfg["E_range"][1], 3000)
    E = np.asarray(E, float)

    norm = n_nucleons(mass_kt) * pot
    rates = {}
    for initial in ("numu", "numubar", "nue", "nuebar"):
        phi = flux(initial, E)
        anti = FLAVOUR_INFO[initial][1]
        sigma = cross_section_cc(E, antineutrino=anti)
        for final in ("e", "mu"):
            key = f"{initial}->{final_flavour(initial, final)}"
            P = probability(osc, initial, final, E, cfg["baseline"], cfg["rho"], matter)
            eff = efficiency.get(key, 1.0) if isinstance(efficiency, dict) else efficiency
            rates[key] = phi * P * sigma * norm * eff
    return E, rates


def group_rates(rates, mode="FHC"):
    """
    Group the 8 channels into the samples an experiment actually looks at.

    FHC: signal = numu->nue,          wrong-sign = numubar->nuebar
    RHC: signal = numubar->nuebar,    wrong-sign = numu->nue
    """
    if mode == "FHC":
        mu, mu_ws, sig, sig_ws = "numu->numu", "numubar->numubar", "numu->nue", "numubar->nuebar"
    else:
        mu, mu_ws, sig, sig_ws = "numubar->numubar", "numu->numu", "numubar->nuebar", "numu->nue"
    return {
        "mu_like":            rates[mu] + rates[mu_ws],
        "mu_like_rightsign":  rates[mu],
        "mu_like_wrongsign":  rates[mu_ws],
        "e_signal":           rates[sig],
        "e_signal_wrongsign": rates[sig_ws],
        "e_beam":             rates["nue->nue"] + rates["nuebar->nuebar"],
        "e_like":             rates[sig] + rates[sig_ws] + rates["nue->nue"] + rates["nuebar->nuebar"],
    }


def integrate(E, dNdE, E_min=-np.inf, E_max=np.inf):
    """Total number of events between E_min and E_max (trapezoidal rule)."""
    m = (E >= E_min) & (E <= E_max)
    return np.trapezoid(dNdE[m], E[m]) if hasattr(np, "trapezoid") else np.trapz(dNdE[m], E[m])


def histogram(E, dNdE, edges):
    """Events per bin for the given bin edges (fine-grid integration)."""
    return np.array([integrate(E, dNdE, lo, hi) for lo, hi in zip(edges[:-1], edges[1:])])


def smear(E, dNdE, E_reco, resolution):
    """
    Toy detector response: Gaussian smearing in energy with sigma = resolution * E.
    Returns dN/dE_reco evaluated on the grid E_reco.
    """
    E = np.asarray(E, float)
    dE = np.gradient(E)
    sigma = resolution * E
    K = np.exp(-0.5 * ((E_reco[:, None] - E[None, :]) / sigma[None, :]) ** 2) \
        / (np.sqrt(2 * np.pi) * sigma[None, :])
    return K @ (dNdE * dE)
