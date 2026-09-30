"""
Default oscillation parameters
==============================

Best-fit values from NuFIT 6.0 (2024), IC24 with SK atmospheric data,
arXiv:2410.05380, www.nu-fit.org.

For each mass ordering we store
    s2_12, s2_13, s2_23 : sin^2 of the mixing angles
    delta_cp            : Dirac CP phase [degrees]
    dm2_21              : Delta m^2_21 [eV^2]  (always > 0)
    dm2_31              : Delta m^2_31 [eV^2]  (> 0 for NO, < 0 for IO)

NuFIT quotes Delta m^2_3l, i.e. Delta m^2_31 for NO and Delta m^2_32 for IO;
we always convert to Delta m^2_31 = Delta m^2_32 + Delta m^2_21.
Feel free to change these numbers: they are just a starting point.
"""

DM2_21 = 7.49e-5

NUFIT_NO = dict(
    s2_12=0.308,
    s2_13=0.02215,
    s2_23=0.470,
    delta_cp=212.0,
    dm2_21=DM2_21,
    dm2_31=+2.513e-3,
)

NUFIT_IO = dict(
    s2_12=0.308,
    s2_13=0.02231,
    s2_23=0.550,
    delta_cp=274.0,
    dm2_21=DM2_21,
    dm2_31=-2.484e-3 + DM2_21,   # Delta m^2_32 + Delta m^2_21
)

DEFAULT_PARAMS = {"NO": NUFIT_NO, "IO": NUFIT_IO}


def default_params(ordering="NO"):
    """Return a *copy* of the default parameters for the given ordering."""
    ordering = ordering.upper()
    if ordering not in DEFAULT_PARAMS:
        raise ValueError("ordering must be 'NO' (normal) or 'IO' (inverted)")
    return dict(DEFAULT_PARAMS[ordering])
