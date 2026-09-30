"""
nuosc - three-flavour neutrino oscillations for students
========================================================

    from nuosc import NeutrinoOscillator
    osc = NeutrinoOscillator(ordering="NO", baseline=1300, rho=2.848)
    osc.matter("mu", "e", E=2.5)      # P(nu_mu -> nu_e) at 2.5 GeV

Units: E [GeV], L [km], dm2 [eV^2], rho [g/cm^3], delta_cp [deg].
"""

from .constants import K_PHASE, K_MATTER, FLAVOURS, FLAVOUR_LABELS
from .parameters import NUFIT_NO, NUFIT_IO, default_params
from .formulas import (pmns_matrix, matter_potential, prob_vacuum,
                       prob_matter_nufast, prob_matter_exact)
from .experiments import EXPERIMENTS, ACCELERATORS, get_experiment
from .fluxes import Flux, load_flux, FLAVOUR_LATEX
from .rates import (expected_events, group_rates, integrate, histogram, smear,
                    cross_section_cc)
from .oscillator import NeutrinoOscillator
from .spectra import load_fd_spectra, predict

__all__ = [
    "NeutrinoOscillator",
    "pmns_matrix", "matter_potential",
    "prob_vacuum", "prob_matter_nufast", "prob_matter_exact",
    "EXPERIMENTS", "ACCELERATORS", "get_experiment",
    "Flux", "load_flux", "FLAVOUR_LATEX",
    "expected_events", "group_rates", "integrate", "histogram", "smear", "cross_section_cc",
    "load_fd_spectra", "predict",
    "NUFIT_NO", "NUFIT_IO", "default_params",
    "K_PHASE", "K_MATTER", "FLAVOURS", "FLAVOUR_LABELS",
]
