"""
Far-detector spectra ready to be multiplied by oscillation probabilities
=========================================================================

The files far_detector_spectra/<EXP>_FD_spectra.csv (EXP = T2K, NOvA, DUNE, HyperK)
contain, per energy bin and beam mode, the number of selected events you would
get if the oscillation probability were 1. The prediction is simply

    mu_like = P(numu->numu) * numu_mu + P(numubar->numubar) * numubar_mu + NC_mu
    e_like  = P(numu->nue)  * numu_e  + P(numubar->nuebar)  * numubar_e
            + P(nue->nue)   * nue_e   + P(nuebar->nuebar)   * nuebar_e
            + P(numu->numu) * numu_misid_e + P(numubar->numubar) * numubar_misid_e + NC_e

>>> from nuosc.spectra import predict
>>> s = predict("DUNE", NeutrinoOscillator("NO"))
>>> s["FHC"]["e_like"].sum()

See far_detector_spectra/README.md for how each file was made.
"""

from pathlib import Path

from .paths import ROOT

import numpy as np

from .formulas import prob_matter_nufast, prob_vacuum

SPECTRA_DIR = ROOT / "far_detector_spectra"

# baseline [km] and density [g/cm^3] used for each experiment
FD_SETUP = {
    "T2K":    dict(baseline=295.0, rho=2.6),
    "HyperK": dict(baseline=295.0, rho=2.6),
    "NOvA":   dict(baseline=810.0, rho=2.84),
    "DUNE":   dict(baseline=1284.9, rho=2.848),
}
COMPONENTS = ("numu_mu", "numubar_mu", "numu_e", "numubar_e", "nue_e", "nuebar_e",
              "numu_misid_e", "numubar_misid_e", "NC_mu", "NC_e")


def load_fd_spectra(name, directory=None):
    """
    Read <name>_FD_spectra.csv. Returns a dict with 'E_lo', 'E_hi', 'E' (bin centres),
    'edges', and for mode in ('FHC', 'RHC') a dict of the components above.
    """
    directory = Path(directory) if directory is not None else SPECTRA_DIR
    path = directory / f"{name}_FD_spectra.csv"
    header = [l[2:].strip() for l in open(path) if l.startswith("#")]
    cols = header[-1].split(",")
    d = np.loadtxt(path, delimiter=",")
    col = {c: d[:, i] for i, c in enumerate(cols)}
    out = dict(name=name, E_lo=col["E_lo"], E_hi=col["E_hi"], E=col["E_center"],
               edges=np.append(col["E_lo"], col["E_hi"][-1]), header="\n".join(header[:-1]),
               **FD_SETUP[name])
    for m in ("FHC", "RHC"):
        out[m] = {c: col[f"{m}_{c}"] for c in COMPONENTS}
    return out


def predict(name, osc, matter=True, spectra=None):
    """
    Oscillated selected spectra for experiment `name` with the parameters of `osc`
    (a NeutrinoOscillator; its baseline/density are NOT used, those of the file are).

    Returns {'FHC': {...}, 'RHC': {...}, 'E': ..., 'edges': ...}; each mode contains
    'mu_like', 'e_like' and every oscillated component (e.g. 'numu_e' = P * numu_e).
    """
    s = load_fd_spectra(name) if spectra is None else spectra
    E, L, rho = s["E"], s["baseline"], s["rho"]
    P = {}
    for anti in (False, True):
        if matter:
            P[anti] = prob_matter_nufast(E, L, **osc.params, rho=rho, Ye=osc.Ye, antineutrino=anti)
        else:
            P[anti] = prob_vacuum(E, L, **osc.params, antineutrino=anti)
    weight = {"numu_mu": P[False][1, 1], "numubar_mu": P[True][1, 1],
              "numu_e": P[False][1, 0], "numubar_e": P[True][1, 0],
              "nue_e": P[False][0, 0], "nuebar_e": P[True][0, 0],
              "numu_misid_e": P[False][1, 1], "numubar_misid_e": P[True][1, 1],
              "NC_mu": 1.0, "NC_e": 1.0}
    out = dict(E=E, edges=s["edges"])
    for m in ("FHC", "RHC"):
        comp = {c: weight[c] * s[m][c] for c in COMPONENTS}
        comp["mu_like"] = comp["numu_mu"] + comp["numubar_mu"] + comp["NC_mu"]
        comp["e_like"] = sum(comp[c] for c in COMPONENTS if c.endswith("_e"))
        out[m] = comp
    return out
