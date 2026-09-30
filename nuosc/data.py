"""
Public far-detector data
========================

NOvA 2024 data release (FERMILAB-DATA-2025-05,
https://publicdocs.fnal.gov/cgi-bin/ShowDocument?docid=594):
26.61e20 POT neutrino beam + 12.5e20 POT antineutrino beam.

The original ROOT files are in data/NOvA_2024_FD_release/data_release_2024/.
They were converted to plain CSV (data/NOvA_2024_FD_release/csv/) so that no
ROOT/uproot installation is needed:

  NOvA_2024_FD_numu_{FHC,RHC}.csv : numu CC sample vs reconstructed energy
      columns E_lo, E_hi, data, osc_total, osc_signal, osc_beam_bkg,
              noosc_total, noosc_signal, noosc_beam_bkg, cosmic_bkg
      (predictions at the NOvA best fit and without oscillations)
  NOvA_2024_FD_nue_totals.csv     : total events of the nue samples

The NOvA "no oscillation" prediction already contains the true FD flux,
cross sections, efficiencies and energy resolution. Multiplying it by our
oscillation probability is therefore the simplest *realistic* prediction we
can make for NOvA (see `reweight_nova_numu`).
"""

from pathlib import Path

from .paths import ROOT

import numpy as np

from .experiments import get_experiment
from .formulas import prob_vacuum, prob_matter_nufast

DATA_DIR = ROOT / "data"
NOVA_2024_CSV = DATA_DIR / "NOvA_2024_FD_release" / "csv"
NOVA_2024_POT = {"FHC": 26.61e20, "RHC": 12.5e20}

_NUMU_COLS = ("E_lo", "E_hi", "data", "osc_total", "osc_signal", "osc_beam_bkg",
              "noosc_total", "noosc_signal", "noosc_beam_bkg", "cosmic_bkg")


def load_nova_2024_numu(mode="FHC", directory=None):
    """NOvA 2024 FD numu CC sample: dict of arrays (see module docstring) + 'edges'."""
    directory = Path(directory) if directory is not None else NOVA_2024_CSV
    d = np.loadtxt(directory / f"NOvA_2024_FD_numu_{mode}.csv", delimiter=",")
    out = {k: d[:, i] for i, k in enumerate(_NUMU_COLS)}
    out["edges"] = np.append(out["E_lo"], out["E_hi"][-1])
    out["pot"] = NOVA_2024_POT[mode]
    return out


def load_nova_2024_nue_totals(directory=None):
    """NOvA 2024 FD nue samples: dict sample -> dict of totals."""
    directory = Path(directory) if directory is not None else NOVA_2024_CSV
    rows = {}
    with open(directory / "NOvA_2024_FD_nue_totals.csv") as f:
        header = None
        for line in f:
            if line.startswith("# sample"):
                header = line[2:].strip().split(",")
            elif not line.startswith("#") and line.strip():
                v = line.strip().split(",")
                rows[v[0]] = {k: float(x) for k, x in zip(header[1:], v[1:])}
    return rows


def smeared_bin_probability(osc, edges, resolution, alpha=1, beta=1,
                            antineutrino=False, matter=True, experiment="NOvA",
                            weights=None, n_true=4000):
    """
    Oscillation probability averaged over reconstructed-energy bins.

    For each reco bin [E_lo, E_hi] we average P(E_true) over the true energies
    that end up in that bin, assuming a Gaussian resolution sigma = resolution*E:

        <P>_bin = int dE w(E) P(E) R_bin(E) / int dE w(E) R_bin(E),
        R_bin(E) = probability that an event of true energy E is reconstructed in the bin,
        w(E)     = unoscillated true spectrum; `weights` = (E_centres, dN/dE) or None (flat).
    """
    from math import erf
    cfg = get_experiment(experiment)
    E = np.linspace(0.05, 3 * edges[-1], n_true)
    if matter:
        P = prob_matter_nufast(E, cfg["baseline"], **osc.params, rho=cfg["rho"],
                               Ye=osc.Ye, antineutrino=antineutrino)[alpha, beta]
    else:
        P = prob_vacuum(E, cfg["baseline"], **osc.params, antineutrino=antineutrino)[alpha, beta]
    w = np.ones_like(E) if weights is None else np.interp(E, *weights, left=0, right=0)
    sigma = resolution * E
    verf = np.vectorize(erf)
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        R = 0.5 * (verf((hi - E) / (np.sqrt(2) * sigma)) - verf((lo - E) / (np.sqrt(2) * sigma)))
        out.append(np.sum(w * P * R) / np.sum(w * R))
    return np.array(out)


def reweight_nova_numu(osc, mode="FHC", resolution=0.09, matter=True, sample=None):
    """
    Our prediction for the NOvA 2024 FD numu sample:

        N_i = N_i^(no osc, signal) * <P_mumu>_i  +  beam background  +  cosmics

    using the oscillation parameters of `osc`. The beam background (mostly NC)
    is taken from NOvA's oscillated prediction.
    """
    s = load_nova_2024_numu(mode) if sample is None else sample
    # use the unoscillated reco spectrum as a proxy for the true spectrum shape
    centres = 0.5 * (s["E_lo"] + s["E_hi"])
    w = (centres, s["noosc_signal"] / (s["E_hi"] - s["E_lo"]))
    P = smeared_bin_probability(osc, s["edges"], resolution, 1, 1,
                                antineutrino=(mode == "RHC"), matter=matter, weights=w)
    return s["noosc_signal"] * P + s["osc_beam_bkg"] + s["cosmic_bkg"]


def poisson_chi2(data, pred):
    """Poisson likelihood-ratio chi^2 = 2 sum [mu - n + n ln(n/mu)]."""
    data = np.asarray(data, float)
    pred = np.asarray(pred, float)
    term = np.where(data > 0, data * np.log(np.where(data > 0, data, 1) / pred), 0.0)
    return 2.0 * np.sum(pred - data + term)


# --------------------------------------------------------------------------- #
#  T2K Run 1-9 (arXiv:2101.03779, PRD 103, 112008)                            #
# --------------------------------------------------------------------------- #
# T2K does not publish its far-detector spectra; the public information used here is
#   * Table XI of the paper: predicted total events per Super-K sample for four values
#     of delta_CP, plus the observed numbers of events;
#   * the data release (data/T2K_OA2019_release/Analysis_*.root, from
#     https://t2k-experiment.org/results/t2kdata-oa-2019/): Delta chi^2 vs delta_CP.
T2K_OA2019_CSV = DATA_DIR / "T2K_OA2019_release" / "csv"
T2K_OA2019_POT = {"FHC": 14.94e20, "RHC": 16.35e20}

# Table III of the paper (reference oscillation parameters)
T2K_OA2019_PARAMS = dict(s2_12=0.304, s2_13=0.0212, s2_23=0.528, dm2_21=7.53e-5,
                         dm2_32_NO=2.509e-3, dm2_31_IO=-2.509e-3, delta_cp_rad=-1.601)

T2K_SAMPLES = ("FHC_mu_like", "RHC_mu_like", "FHC_e_like", "RHC_e_like", "FHC_e_CC1pi")


def t2k_oscillator(delta_cp_rad=T2K_OA2019_PARAMS["delta_cp_rad"], ordering="NO"):
    """NeutrinoOscillator with the T2K Table III parameters (delta_CP in radians)."""
    from .oscillator import NeutrinoOscillator
    p = T2K_OA2019_PARAMS
    dm2_31 = p["dm2_32_NO"] + p["dm2_21"] if ordering == "NO" else p["dm2_31_IO"]
    return NeutrinoOscillator(ordering, s2_12=p["s2_12"], s2_13=p["s2_13"], s2_23=p["s2_23"],
                              dm2_21=p["dm2_21"], dm2_31=dm2_31,
                              delta_cp=np.degrees(delta_cp_rad))


def load_t2k_oa2019_rates(directory=None):
    """
    T2K Table XI. Returns (predictions, observed):
      predictions : dict delta_cp_rad -> dict sample -> events
      observed    : dict sample -> events
    """
    directory = Path(directory) if directory is not None else T2K_OA2019_CSV
    pred, obs = {}, {}
    for line in open(directory / "T2K_OA2019_SK_event_rates.csv"):
        if line.startswith("#") or not line.strip():
            continue
        v = line.strip().split(",")
        vals = dict(zip(T2K_SAMPLES, map(float, v[1:])))
        if v[0] == "observed":
            obs = vals
        else:
            pred[float(v[0])] = vals
    return pred, obs


def load_t2k_oa2019_dchi2_dcp(ordering="NO", directory=None):
    """Official T2K (analysis A, with reactor constraint) Delta chi^2 vs delta_CP [rad]."""
    directory = Path(directory) if directory is not None else T2K_OA2019_CSV
    d = np.loadtxt(directory / f"T2K_OA2019_A_dchi2_dcp_wRC_{ordering}.csv", delimiter=",")
    return d[:, 0], d[:, 1]
