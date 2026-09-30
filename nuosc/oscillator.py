"""
NeutrinoOscillator: a small user-friendly wrapper
=================================================

Example
-------
>>> from nuosc import NeutrinoOscillator
>>> osc = NeutrinoOscillator(ordering="NO", baseline=1300, rho=2.848)
>>> E = np.linspace(0.5, 5, 500)                  # GeV
>>> P_mue_vac = osc.vacuum("mu", "e", E)          # P(nu_mu -> nu_e) in vacuum
>>> P_mue_mat = osc.matter("mu", "e", E)          # ... in matter (NuFast)
>>> P_all     = osc.matter(E=E)                   # full 3x3 matrix, shape (3,3,500)
>>> osc.ordering = "IO"                           # switch to inverted ordering
>>> osc.set_params(delta_cp=270, s2_23=0.5)       # change any parameter
"""

import numpy as np

from .constants import FLAVOURS
from .parameters import default_params
from .experiments import get_experiment
from .formulas import (pmns_matrix, prob_vacuum, prob_matter_nufast,
                       prob_matter_exact)

PARAM_NAMES = ("s2_12", "s2_13", "s2_23", "delta_cp", "dm2_21", "dm2_31")


def _flavour_index(f):
    """'e'/'mu'/'tau' (or 0/1/2) -> 0/1/2."""
    if isinstance(f, (int, np.integer)):
        if f in (0, 1, 2):
            return int(f)
    elif isinstance(f, str) and f.lower() in FLAVOURS:
        return FLAVOURS[f.lower()]
    raise ValueError(f"Unknown flavour {f!r}: use 'e', 'mu', 'tau' or 0, 1, 2")


class NeutrinoOscillator:
    """
    Three-flavour neutrino oscillations in vacuum and in constant-density matter.

    Parameters
    ----------
    ordering     : 'NO' (normal) or 'IO' (inverted); loads NuFIT 6.0 best fit
    baseline     : L [km]
    energy       : default energy grid E [GeV] (scalar or array)
    rho          : constant matter density [g/cm^3]
    Ye           : electron fraction (0.5 for the Earth's crust/mantle)
    antineutrino : False for neutrinos, True for antineutrinos
    N_Newton     : Newton iterations in NuFast (0 is already very accurate)
    **params     : override any of s2_12, s2_13, s2_23, delta_cp [deg],
                   dm2_21 [eV^2], dm2_31 [eV^2]
    """

    def __init__(self, ordering="NO", baseline=1300.0, energy=None, rho=2.8,
                 Ye=0.5, antineutrino=False, N_Newton=0, **params):
        self.baseline = baseline
        self.energy = np.linspace(0.2, 10.0, 2000) if energy is None else energy
        self.rho = rho
        self.Ye = Ye
        self.antineutrino = antineutrino
        self.N_Newton = N_Newton
        self.experiment = None            # set by from_experiment()

        self.ordering = ordering          # loads default parameters
        self.set_params(**params)

    @classmethod
    def from_experiment(cls, name, ordering="NO", n_energy=2000, **kwargs):
        """
        Build an oscillator from an entry of nuosc.experiments.EXPERIMENTS,
        e.g. NeutrinoOscillator.from_experiment("NOvA", ordering="IO").

        Any keyword (antineutrino=True, delta_cp=0, rho=3.0, ...) overrides
        the experiment defaults. The experiment config is kept in self.experiment.
        """
        cfg = get_experiment(name)
        setup = dict(baseline=cfg["baseline"], rho=cfg["rho"],
                     antineutrino=cfg["antineutrino"],
                     energy=np.linspace(*cfg["E_range"], n_energy))
        setup.update(kwargs)
        osc = cls(ordering=ordering, **setup)
        osc.experiment = cfg
        return osc

    # ------------------------------------------------------------------ #
    #  Parameters                                                        #
    # ------------------------------------------------------------------ #
    @property
    def ordering(self):
        """Mass ordering, 'NO' or 'IO'. Setting it reloads the defaults."""
        return self._ordering

    @ordering.setter
    def ordering(self, value):
        value = value.upper()
        for name, val in default_params(value).items():
            setattr(self, name, val)
        self._ordering = value

    def set_params(self, **params):
        """Change one or more oscillation parameters, e.g. set_params(delta_cp=0)."""
        for name, val in params.items():
            if name == "dm2_32":                     # convenience
                self.dm2_31 = val + self.dm2_21
                continue
            if name not in PARAM_NAMES:
                raise KeyError(f"Unknown parameter {name!r}. "
                               f"Valid: {PARAM_NAMES + ('dm2_32',)}")
            setattr(self, name, float(val))
        if np.sign(self.dm2_31) != (1 if self._ordering == "NO" else -1):
            print(f"Warning: dm2_31 = {self.dm2_31:.3e} has the 'wrong' sign "
                  f"for ordering {self._ordering}.")
        return self

    def reset(self):
        """Go back to the default (NuFIT) values for the current ordering."""
        self.ordering = self._ordering
        return self

    @property
    def params(self):
        """Current oscillation parameters as a dict."""
        return {name: getattr(self, name) for name in PARAM_NAMES}

    @property
    def dm2_32(self):
        return self.dm2_31 - self.dm2_21

    @property
    def theta_12(self):
        return np.rad2deg(np.arcsin(np.sqrt(self.s2_12)))

    @property
    def theta_13(self):
        return np.rad2deg(np.arcsin(np.sqrt(self.s2_13)))

    @property
    def theta_23(self):
        return np.rad2deg(np.arcsin(np.sqrt(self.s2_23)))

    @property
    def U(self):
        """PMNS matrix (3x3 complex)."""
        return pmns_matrix(self.s2_12, self.s2_13, self.s2_23, self.delta_cp)

    def __repr__(self):
        nu = "antineutrino" if self.antineutrino else "neutrino"
        return (f"NeutrinoOscillator({self._ordering}, {nu}, "
                f"L = {self.baseline} km, rho = {self.rho} g/cm^3)\n"
                f"  sin^2(th12) = {self.s2_12:.4f}   (th12 = {self.theta_12:.2f} deg)\n"
                f"  sin^2(th13) = {self.s2_13:.5f}  (th13 = {self.theta_13:.2f} deg)\n"
                f"  sin^2(th23) = {self.s2_23:.4f}   (th23 = {self.theta_23:.2f} deg)\n"
                f"  delta_CP    = {self.delta_cp:.1f} deg\n"
                f"  dm2_21      = {self.dm2_21:.3e} eV^2\n"
                f"  dm2_31      = {self.dm2_31:+.4e} eV^2\n"
                f"  dm2_32      = {self.dm2_32:+.4e} eV^2")

    # ------------------------------------------------------------------ #
    #  Probabilities                                                     #
    # ------------------------------------------------------------------ #
    def probability(self, alpha=None, beta=None, E=None, L=None,
                    matter=True, method="nufast"):
        """
        P(nu_alpha -> nu_beta).

        alpha, beta : 'e', 'mu', 'tau' (or 0, 1, 2). If both are None the
                      full matrix P[alpha, beta, ...] is returned.
        E           : energy [GeV]; default self.energy
        L           : baseline [km]; default self.baseline
        matter      : False -> vacuum formula, True -> constant density rho
        method      : 'nufast' (default) or 'exact' (numerical diagonalisation)
        """
        E = self.energy if E is None else E
        L = self.baseline if L is None else L

        if not matter:
            P = prob_vacuum(E, L, **self.params, antineutrino=self.antineutrino)
        elif method == "nufast":
            P = prob_matter_nufast(E, L, **self.params, rho=self.rho, Ye=self.Ye,
                                   antineutrino=self.antineutrino,
                                   N_Newton=self.N_Newton)
        elif method == "exact":
            P = prob_matter_exact(E, L, **self.params, rho=self.rho, Ye=self.Ye,
                                  antineutrino=self.antineutrino)
        else:
            raise ValueError("method must be 'nufast' or 'exact'")

        if alpha is None and beta is None:
            return P
        return P[_flavour_index(alpha), _flavour_index(beta)]

    def vacuum(self, alpha=None, beta=None, E=None, L=None):
        """Vacuum probability (shortcut for probability(..., matter=False))."""
        return self.probability(alpha, beta, E, L, matter=False)

    def matter(self, alpha=None, beta=None, E=None, L=None, method="nufast"):
        """Constant-density matter probability (NuFast by default)."""
        return self.probability(alpha, beta, E, L, matter=True, method=method)
