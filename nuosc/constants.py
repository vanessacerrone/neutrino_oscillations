"""
Physical constants and unit conventions
=======================================

Units used everywhere in this package
-------------------------------------
    energy E          : GeV      (reactor neutrinos: 3 MeV = 3e-3 GeV)
    baseline L        : km
    mass splittings   : eV^2
    matter density    : g/cm^3
    CP phase delta_cp : degrees  (converted to radians internally)
"""

import numpy as np

# Fermi constant [GeV^-2]  (PDG)
G_F = 1.1663787e-5

# hbar * c [GeV cm]
HBARC_GEV_CM = 1.973269804e-14

# hbar * c [eV km]  ->  1 eV^-1 = 1.97327e-10 km
HBARC_EV_KM = 1.973269804e-10

# Avogadro number [1/mol]  (number of nucleons per gram of matter)
N_A = 6.02214076e23

# ---------------------------------------------------------------------------
# Oscillation phase:   Delta_ij = dm2_ij L / (4 E)   in natural units.
# With dm2 in eV^2, L in km and E in GeV:
#       Delta_ij = K_PHASE * dm2_ij [eV^2] * L [km] / E [GeV]
# K_PHASE = 1e-9 / (4 * hbarc[eV km])  ~ 1.26693   (the famous "1.27")
# ---------------------------------------------------------------------------
K_PHASE = 1e-9 / (4.0 * HBARC_EV_KM)

# ---------------------------------------------------------------------------
# Matter potential (Wolfenstein):   V_CC = sqrt(2) G_F N_e,
# and in the Hamiltonian it appears as  A = 2 E V_CC   [eV^2].
# With N_e = Y_e * rho * N_A (electrons per cm^3):
#       A [eV^2] = K_MATTER * Y_e * rho [g/cm^3] * E [GeV]
# K_MATTER ~ 1.5259e-4
# ---------------------------------------------------------------------------
K_MATTER = 2.0 * np.sqrt(2.0) * G_F * HBARC_GEV_CM**3 * N_A * 1e18

# Flavour labels <-> indices
FLAVOURS = {"e": 0, "mu": 1, "tau": 2}
FLAVOUR_LABELS = {0: r"e", 1: r"\mu", 2: r"\tau"}
