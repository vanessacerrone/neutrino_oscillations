"""
Build the easy-to-use far-detector spectra in far_detector_spectra/.

    python tools/build_fd_spectra.py

For every experiment and beam mode (FHC = neutrino beam, RHC = antineutrino beam)
the output CSV gives, in bins of TRUE neutrino energy, the number of SELECTED
events you would get if the oscillation probability were 1:

    <mode>_numu_mu     : nu_mu flux    -> nu_mu CC,    mu-like sample   (x P(numu->numu))
    <mode>_numubar_mu  : nu_mu-bar flux-> nu_mu-bar CC, mu-like sample  (x P(numubar->numubar))
    <mode>_numu_e      : nu_mu flux    -> nu_e CC,     e-like sample    (x P(numu->nue))
    <mode>_numubar_e   : nu_mu-bar flux-> nu_e-bar CC, e-like sample    (x P(numubar->nuebar))
    <mode>_nue_e       : intrinsic nu_e     -> e-like sample            (x P(nue->nue))
    <mode>_nuebar_e    : intrinsic nu_e-bar -> e-like sample            (x P(nuebar->nuebar))
    <mode>_numu_misid_e, <mode>_numubar_misid_e : nu_mu CC mis-identified as e-like (x P(numu->numu));
                         only available for DUNE (zero elsewhere)
    <mode>_NC_mu, <mode>_NC_e : neutral currents (+ cosmics for NOvA), NOT oscillated

Ingredients (best public information, see far_detector_spectra/README.md):
  DUNE   : official TDR GLoBES configuration (flux, GENIE xsec, smearing, efficiencies),
           selected events projected back on true energy.
  HyperK : T2K/SK flux (same beam, same 2.5 deg off-axis angle and 295 km), 187 kt,
           GENIE sigma/E (smoothed), constant selection efficiencies per component
           calibrated to Tables XXXVI-XXXVII of the Hyper-K Design Report (arXiv:1805.04163).
  T2K    : same as Hyper-K (Super-K = same detector technology) scaled to 22.5 kt and the
           Run 1-9 exposure, then one normalisation per sample adjusted to Table XI of
           arXiv:2101.03779.
  NOvA   : mu-like: NOvA's own no-oscillation FD prediction (2024 data release), split into
           nu / nu-bar with our flux model; e-like: same FD nu_mu spectrum shape, normalised
           to Table 2 of arXiv:2509.04361 at the NOvA best fit.
"""

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from nuosc import NeutrinoOscillator, load_flux                       # noqa: E402
from nuosc.formulas import prob_matter_nufast                          # noqa: E402
from nuosc.dune_globes import DUNEGlobes, DUNE_GLOBES_DIR              # noqa: E402
from nuosc.data import load_nova_2024_numu                             # noqa: E402

OUT = ROOT / "far_detector_spectra"
N_A = 6.02214076e23
COLS = ["numu_mu", "numubar_mu", "numu_e", "numubar_e", "nue_e", "nuebar_e",
        "numu_misid_e", "numubar_misid_e", "NC_mu", "NC_e"]


# --------------------------------------------------------------------------- #
#  helpers                                                                    #
# --------------------------------------------------------------------------- #
def smoothed_genie():
    """GENIE sigma/E tables (1e-38 cm^2/GeV/nucleon) from the DUNE config, smoothed in log E."""
    out = {}
    for kind in ("cc", "nc"):
        d = np.loadtxt(DUNE_GLOBES_DIR / "xsec" / f"xsec_{kind}.dat")
        k = np.ones(7) / 7
        sm = d.copy()
        for j in range(1, 7):
            sm[:, j] = np.convolve(np.pad(d[:, j], 3, mode="edge"), k, mode="valid")
        out[kind] = sm
    return out


XS = smoothed_genie()
_XCOL = {"nue": 1, "numu": 2, "nuebar": 4, "numubar": 5}


def sigma(kind, flav, E):
    """cross section per nucleon [cm^2]"""
    d = XS[kind]
    return np.interp(np.log10(np.maximum(E, 1e-3)), d[:, 0], d[:, _XCOL[flav]], left=0, right=d[-1, _XCOL[flav]]) * E * 1e-38


def binned_rate(flux, flav_flux, flav_xs, kind, edges, mass_kt, pot, n_sub=20):
    """Phi * sigma * N * POT integrated in each bin (no oscillation, no efficiency)."""
    out = np.zeros(len(edges) - 1)
    for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        e = np.linspace(lo, hi, n_sub + 1)
        e = 0.5 * (e[1:] + e[:-1])
        out[i] = np.sum(flux(flav_flux, e) * sigma(kind, flav_xs, e)) * (hi - lo) / n_sub
    return out * mass_kt * 1e9 * N_A * pot


def probs(osc, E, L, rho):
    """P[alpha, beta] for neutrinos and antineutrinos at energies E (matter)."""
    return {anti: prob_matter_nufast(E, L, **osc.params, rho=rho, antineutrino=anti)
            for anti in (False, True)}


def combine(table, mode, osc, E, L, rho):
    """mu-like and e-like totals from the table columns (the recipe in the README)."""
    P = probs(osc, E, L, rho)
    t = {c: table[f"{mode}_{c}"] for c in COLS}
    mu = P[False][1, 1] * t["numu_mu"] + P[True][1, 1] * t["numubar_mu"] + t["NC_mu"]
    e = (P[False][1, 0] * t["numu_e"] + P[True][1, 0] * t["numubar_e"]
         + P[False][0, 0] * t["nue_e"] + P[True][0, 0] * t["nuebar_e"]
         + P[False][1, 1] * t["numu_misid_e"] + P[True][1, 1] * t["numubar_misid_e"] + t["NC_e"])
    return mu, e


def write(name, edges, table, header):
    E = 0.5 * (edges[1:] + edges[:-1])
    cols = ["E_lo", "E_hi", "E_center"] + [f"{m}_{c}" for m in ("FHC", "RHC") for c in COLS]
    data = np.column_stack([edges[:-1], edges[1:], E] + [table[c] for c in cols[3:]])
    OUT.mkdir(exist_ok=True)
    np.savetxt(OUT / f"{name}_FD_spectra.csv", data, delimiter=",", fmt="%.6g",
               header=header.strip() + "\n" + ",".join(cols), comments="# ")
    print(f"wrote {name}_FD_spectra.csv ({len(E)} bins)")


# --------------------------------------------------------------------------- #
#  DUNE                                                                       #
# --------------------------------------------------------------------------- #
def build_dune():
    d = DUNEGlobes()
    osc0 = None
    edges = d.sampling_edges                      # true energy 0-110 GeV (GLoBES sampling bins)
    n = len(edges) - 1
    table = {}

    def proj(ch):
        c = d.channels[ch]
        n_true = d.channel_true(osc0, ch)
        if c["oscillate"] and c["initial"] != c["final"]:
            # appearance channel: channel_true(None) is zero -> recompute with P = 1
            E = d.E_true
            n_true = (d.normalisation(c["mode"]) * d.flux_at(c["mode"], c["initial"], c["polarity"], E)
                      * d.xsec_at(c["xsec"], c["final"], c["polarity"], E) * d.dE_true)
        rule = [r for r, v in d.rules.items() if ch in v["signal"] + v["background"]][0]
        mask = d.sample(None, rule)["mask"]
        sel = (c["smear"] * (c["eff"] * mask)[:, None]).sum(axis=0) * n_true
        return sel[:n]

    for m in ("FHC", "RHC"):
        table[f"{m}_numu_mu"] = proj(f"{m}_dis_sig_numu")
        table[f"{m}_numubar_mu"] = proj(f"{m}_dis_sig_numubar")
        table[f"{m}_NC_mu"] = proj(f"{m}_dis_bkg_nuNC") + proj(f"{m}_dis_bkg_nubarNC")
        table[f"{m}_numu_e"] = proj(f"{m}_app_osc_nue")
        table[f"{m}_numubar_e"] = proj(f"{m}_app_osc_nuebar")
        table[f"{m}_nue_e"] = proj(f"{m}_app_bkg_nue")
        table[f"{m}_nuebar_e"] = proj(f"{m}_app_bkg_nuebar")
        table[f"{m}_numu_misid_e"] = proj(f"{m}_app_bkg_numu")
        table[f"{m}_numubar_misid_e"] = proj(f"{m}_app_bkg_numubar")
        table[f"{m}_NC_e"] = proj(f"{m}_app_bkg_nuNC") + proj(f"{m}_app_bkg_nubarNC")

    # validation against the full GLoBES calculation (reco energy)
    E = 0.5 * (edges[1:] + edges[:-1])
    osc = NeutrinoOscillator("NO")
    S = d.all_samples(osc)
    for m, (app, dis) in {"FHC": ("nue_app", "numu_dis"), "RHC": ("nuebar_app", "numubar_dis")}.items():
        mu, e = combine(table, m, osc, E, d.baseline, d.rho)
        print(f"  DUNE {m}: mu-like {mu.sum():8.1f} (GLoBES {S[dis]['total'].sum():8.1f}),"
              f" e-like {e.sum():7.1f} (GLoBES {S[app]['total'].sum():7.1f})")

    header = f"""
DUNE far detector: selected events per TRUE-energy bin if P = 1 (see far_detector_spectra/README.md).
Source: DUNE TDR GLoBES configuration, arXiv:2103.04797 (flux, GENIE cross sections on Ar, smearing,
post-smearing efficiencies; selected events projected back on true energy, reco window 0.5-18 GeV).
Exposure: {d.years['FHC']:g} yr FHC + {d.years['RHC']:g} yr RHC x {d.power:g}e20 POT/yr, {d.mass_kt:g} kt (624 kt-MW-yr).
Baseline {d.baseline} km, density {d.rho} g/cm^3. nu_tau CC events are not included (zero cross section in the config).
Recipe: mu_like = P_mumu*numu_mu + Pbar_mumu*numubar_mu + NC_mu ;
        e_like  = P_mue*numu_e + Pbar_mue*numubar_e + P_ee*nue_e + Pbar_ee*nuebar_e
                  + P_mumu*numu_misid_e + Pbar_mumu*numubar_misid_e + NC_e
Energies in GeV, events per bin."""
    write("DUNE", edges, table, header)


# --------------------------------------------------------------------------- #
#  Hyper-K and T2K (Super-K technology)                                       #
# --------------------------------------------------------------------------- #
# Hyper-K Design Report, arXiv:1805.04163, Tables XXXVI and XXXVII (NO, sin^2 2th13 = 0.1,
# delta = 0, sin^2 th23 = 0.5, dm2_32 = 2.4e-3, sin^2 2th12 = 0.8704, dm2_21 = 7.6e-5)
HK_DR = {
    "FHC": dict(numu_e=1643, numubar_e=15, nue_e=248, nuebar_e=11, NC_e=134,
                numu_mu=6043 + 2981, numubar_mu=348 + 194, NC_mu=480),
    "RHC": dict(numu_e=206, numubar_e=1183, nue_e=101, nuebar_e=216, NC_e=196,
                numu_mu=2699 + 2354, numubar_mu=6099 + 1961, NC_mu=603),
}
HK_POT = {"FHC": 2.7e22 / 4, "RHC": 2.7e22 * 3 / 4}      # 2.7e22 POT, nu : nubar = 1 : 3
HK_MASS = 187.0
E_CUT_E = 1.25                                          # e-like: E_rec < 1.25 GeV


def hk_dr_oscillator():
    s2_13 = (1 - np.sqrt(1 - 0.10)) / 2
    s2_12 = (1 - np.sqrt(1 - 0.8704)) / 2
    return NeutrinoOscillator("NO", s2_12=s2_12, s2_13=s2_13, s2_23=0.5, delta_cp=0.0,
                              dm2_21=7.6e-5, dm2_31=2.4e-3 + 7.6e-5)


def sk_like_table(mass_kt, pot, edges):
    """Raw (efficiency 1) SK-type rates with the T2K flux, before calibration."""
    raw = {}
    for m in ("FHC", "RHC"):
        f = load_flux("T2K", m)
        f.interpolate = False
        r = lambda ff, fx, kind="cc": binned_rate(f, ff, fx, kind, edges, mass_kt, pot[m])
        raw[m] = dict(numu_mu=r("numu", "numu"), numubar_mu=r("numubar", "numubar"),
                      numu_e=r("numu", "nue"), numubar_e=r("numubar", "nuebar"),
                      nue_e=r("nue", "nue"), nuebar_e=r("nuebar", "nuebar"),
                      NC=r("numu", "numu", "nc") + r("numubar", "numubar", "nc")
                      + r("nue", "nue", "nc") + r("nuebar", "nuebar", "nc"))
    return raw


def calibrate_sk(raw, edges, targets, osc, L=295.0, rho=2.6):
    """Constant efficiency per component so that the totals match `targets` for `osc`."""
    E = 0.5 * (edges[1:] + edges[:-1])
    P = probs(osc, E, L, rho)
    ecut = (E < E_CUT_E).astype(float)
    table, eff = {}, {}
    for m in ("FHC", "RHC"):
        pw = {"numu_mu": P[False][1, 1], "numubar_mu": P[True][1, 1], "numu_e": P[False][1, 0],
              "numubar_e": P[True][1, 0], "nue_e": P[False][0, 0], "nuebar_e": P[True][0, 0]}
        for c, w in pw.items():
            shape = raw[m][c] * (ecut if c.endswith("_e") else 1.0)
            eff[(m, c)] = targets[m][c] / np.sum(w * shape)
            table[f"{m}_{c}"] = eff[(m, c)] * shape
        table[f"{m}_numu_misid_e"] = np.zeros(len(E))
        table[f"{m}_numubar_misid_e"] = np.zeros(len(E))
        table[f"{m}_NC_mu"] = targets[m]["NC_mu"] * raw[m]["NC"] / raw[m]["NC"].sum()
        nc_e = raw[m]["NC"] * ecut
        table[f"{m}_NC_e"] = targets[m]["NC_e"] * nc_e / nc_e.sum()
    return table, eff


def build_hyperk_and_t2k():
    f = load_flux("T2K", "FHC")
    edges = f.edges["numu"]
    edges = edges[edges <= 10.0001]
    # --- Hyper-K ---
    raw = sk_like_table(HK_MASS, HK_POT, edges)
    table, eff = calibrate_sk(raw, edges, HK_DR, hk_dr_oscillator())
    print("  HyperK effective efficiencies:", {f"{m}_{c}": round(v, 3) for (m, c), v in eff.items()})
    E = 0.5 * (edges[1:] + edges[:-1])
    for m in ("FHC", "RHC"):
        mu, e = combine(table, m, hk_dr_oscillator(), E, 295.0, 2.6)
        print(f"  HyperK {m} at DR parameters: mu-like {mu.sum():.0f}, e-like {e.sum():.0f}"
              f" (DR: {sum(HK_DR[m][k] for k in ('numu_mu','numubar_mu','NC_mu'))},"
              f" {sum(HK_DR[m][k] for k in ('numu_e','numubar_e','nue_e','nuebar_e','NC_e'))})")
        mu, e = combine(table, m, NeutrinoOscillator("NO"), E, 295.0, 2.6)
        print(f"  HyperK {m} NuFIT 6.0 NO: mu-like {mu.sum():.0f}, e-like {e.sum():.0f}")
    header = f"""
Hyper-Kamiokande far detector: selected events per TRUE-energy bin if P = 1 (see far_detector_spectra/README.md).
Flux: T2K 2020 flux at Super-K (same J-PARC beam, 2.5 deg off-axis, 295 km). Cross section: GENIE sigma/E (smoothed).
Selection: constant efficiency per component, calibrated to Tables XXXVI-XXXVII of the Hyper-K Design Report
(arXiv:1805.04163); e-like columns only below E = {E_CUT_E} GeV (E_rec < 1.25 GeV cut); NC shapes approximate.
Exposure: 2.7e22 POT (1.3 MW x 10 yr), nu:nubar = 1:3 -> FHC {HK_POT['FHC']:.3g} POT, RHC {HK_POT['RHC']:.3g} POT; {HK_MASS:g} kt fiducial (1 tank).
Baseline 295 km, density 2.6 g/cm^3.
Recipe: mu_like = P_mumu*numu_mu + Pbar_mumu*numubar_mu + NC_mu ;
        e_like  = P_mue*numu_e + Pbar_mue*numubar_e + P_ee*nue_e + Pbar_ee*nuebar_e
                  + P_mumu*numu_misid_e + Pbar_mumu*numubar_misid_e + NC_e
Energies in GeV, events per bin."""
    write("HyperK", edges, table, header)

    # --- T2K: same efficiencies, SK mass and Run 1-9 exposure ---
    t2k_pot = {"FHC": 14.94e20, "RHC": 16.35e20}
    scale = {m: (22.5 / HK_MASS) * t2k_pot[m] / HK_POT[m] for m in ("FHC", "RHC")}
    t2k = {k: v * scale[k[:3]] for k, v in table.items()}
    # one normalisation per SK sample to match Table XI of arXiv:2101.03779 at delta = -pi/2
    from nuosc.data import t2k_oscillator, load_t2k_oa2019_rates
    pred, obs = load_t2k_oa2019_rates()
    o = t2k_oscillator(-1.5708)
    norm = {}
    for m, smu, se in (("FHC", "FHC_mu_like", "FHC_e_like"), ("RHC", "RHC_mu_like", "RHC_e_like")):
        mu, e = combine(t2k, m, o, E, 295.0, 2.6)
        norm[(m, "mu")] = pred[-1.5708][smu] / mu.sum()
        norm[(m, "e")] = pred[-1.5708][se] / e.sum()
        for c in COLS:
            t2k[f"{m}_{c}"] = t2k[f"{m}_{c}"] * norm[(m, "mu" if c.endswith("_mu") else "e")]
    print("  T2K sample normalisations (Table XI / HyperK-scaled):", {f"{k[0]}_{k[1]}": round(v, 3) for k, v in norm.items()})
    for d in sorted(pred):
        vals = []
        for m, smu, se in (("FHC", "FHC_mu_like", "FHC_e_like"), ("RHC", "RHC_mu_like", "RHC_e_like")):
            mu, e = combine(t2k, m, t2k_oscillator(d), E, 295.0, 2.6)
            vals += [f"{m} mu {mu.sum():6.1f} ({pred[d][smu]:6.1f})", f"e {e.sum():5.1f} ({pred[d][se]:5.1f})"]
        print(f"  T2K delta = {d:+.2f}: " + ", ".join(vals))
    header = f"""
T2K far detector (Super-K): selected events per TRUE-energy bin if P = 1 (see far_detector_spectra/README.md).
Built like HyperK_FD_spectra.csv (T2K 2020 SK flux, GENIE sigma/E, Hyper-K DR efficiencies) scaled to 22.5 kt and
the Run 1-9 exposure (14.94e20 POT FHC, 16.35e20 POT RHC), then one normalisation per sample (mu-like, e-like; FHC, RHC)
so that the totals match Table XI of arXiv:2101.03779 at delta_CP = -pi/2 (Table III parameters).
The single-ring e-like + 1 decay-electron (CC1pi) sample is not included.
Baseline 295 km, density 2.6 g/cm^3.
Recipe: mu_like = P_mumu*numu_mu + Pbar_mumu*numubar_mu + NC_mu ;
        e_like  = P_mue*numu_e + Pbar_mue*numubar_e + P_ee*nue_e + Pbar_ee*nuebar_e
                  + P_mumu*numu_misid_e + Pbar_mumu*numubar_misid_e + NC_e
Energies in GeV, events per bin."""
    write("T2K", edges, t2k, header)


# --------------------------------------------------------------------------- #
#  NOvA                                                                       #
# --------------------------------------------------------------------------- #
# arXiv:2509.04361 (PRL 136, 011802), Table 2, frequentist best fit with Daya Bay 1D constraint
NOVA_T2 = {
    "FHC": dict(numu_mu=372.3, numubar_mu=24.5, numu_e=125.3, numubar_e=1.8, nue_e=26.1, NC_e=16.8 + 5.5 + 0.8),
    "RHC": dict(numu_mu=24.4, numubar_mu=71.5, numu_e=2.1, numubar_e=18.9, nue_e=6.5, NC_e=2.0 + 1.1 + 0.1),
}


def nova_bf_oscillator():
    s2_13 = (1 - np.sqrt(1 - 0.0851)) / 2          # Daya Bay constraint used by NOvA
    s2_12 = (1 - np.sqrt(1 - 0.851)) / 2           # PDG 2019
    return NeutrinoOscillator("NO", s2_12=s2_12, s2_13=s2_13, s2_23=0.55, delta_cp=0.87 * 180,
                              dm2_21=7.53e-5, dm2_31=2.441e-3 + 7.53e-5)


def build_nova():
    L, rho = 810.0, 2.84
    s = {m: load_nova_2024_numu(m) for m in ("FHC", "RHC")}
    edges = s["FHC"]["edges"]
    E = 0.5 * (edges[1:] + edges[:-1])
    osc = nova_bf_oscillator()
    P = probs(osc, E, L, rho)
    table = {}
    for m in ("FHC", "RHC"):
        f = load_flux("NOvA", m)
        pot = s[m]["pot"]
        r = lambda ff, fx, kind="cc": binned_rate(f, ff, fx, kind, edges, 14.0, pot)
        rnu, rnub = r("numu", "numu"), r("numubar", "numubar")
        frac_nu = rnu / np.maximum(rnu + rnub, 1e-30)
        base = s[m]["noosc_signal"]                    # NOvA FD no-oscillation numu CC (nu + nubar)
        t = NOVA_T2[m]
        # mu-like: NOvA prediction split into nu / nubar, normalised to Table 2
        for c, frac, w in (("numu_mu", frac_nu, P[False][1, 1]), ("numubar_mu", 1 - frac_nu, P[True][1, 1])):
            shape = base * frac
            table[f"{m}_{c}"] = t[c] * shape / np.sum(w * shape)
        table[f"{m}_NC_mu"] = s[m]["osc_beam_bkg"] + s[m]["cosmic_bkg"]
        table[f"{m}_numu_misid_e"] = np.zeros(len(E))
        table[f"{m}_numubar_misid_e"] = np.zeros(len(E))
        # e-like: FD nu_mu CC shape (true FD flux x xsec x efficiency) times a constant
        for c, shape, w in (("numu_e", base * frac_nu, P[False][1, 0]),
                            ("numubar_e", base * (1 - frac_nu), P[True][1, 0])):
            table[f"{m}_{c}"] = t[c] * shape / np.sum(w * shape)
        # beam nu_e: shape from the (ND-scaled) beam nu_e flux, normalised to Table 2
        be = r("nue", "nue") + 1e-30
        beb = r("nuebar", "nuebar")
        tot = np.sum(P[False][0, 0] * be + P[True][0, 0] * beb)
        table[f"{m}_nue_e"] = t["nue_e"] * be / tot
        table[f"{m}_nuebar_e"] = t["nue_e"] * beb / tot
        # NC + cosmics in the e-like sample: shape of the no-oscillation nu_mu spectrum (approximation)
        table[f"{m}_NC_e"] = t["NC_e"] * base / base.sum()
        mu, e = combine(table, m, osc, E, L, rho)
        print(f"  NOvA {m} at NOvA best fit: mu-like {mu.sum():.1f} (NOvA pred {s[m]['osc_total'].sum():.1f},"
              f" data {s[m]['data'].sum():.0f}), e-like {e.sum():.1f}")
        mu, e = combine(table, m, NeutrinoOscillator("NO"), E, L, rho)
        print(f"  NOvA {m} NuFIT 6.0 NO: mu-like {mu.sum():.1f}, e-like {e.sum():.1f}")
    header = """
NOvA far detector: selected events per energy bin if P = 1 (see far_detector_spectra/README.md).
mu-like: NOvA no-oscillation FD prediction of the 2024 data release (reconstructed energy, used as a proxy of
true energy), split into nu / nubar with our flux model; NC_mu = NOvA beam background + cosmics (per bin).
e-like: same FD nu_mu CC spectrum shape (true FD flux x xsec x efficiency), beam nu_e shape from the ND-scaled
flux, NC + cosmics with the nu_mu shape; all normalised to Table 2 of arXiv:2509.04361 (NOvA best fit: NO,
dm2_32 = 2.441e-3, sin^2 th23 = 0.55, delta = 0.87 pi, sin^2 2th13 = 0.0851). Shapes of e-like columns are approximate.
Exposure: 26.61e20 POT FHC, 12.5e20 POT RHC (14 kt). Baseline 810 km, density 2.84 g/cm^3.
Recipe: mu_like = P_mumu*numu_mu + Pbar_mumu*numubar_mu + NC_mu ;
        e_like  = P_mue*numu_e + Pbar_mue*numubar_e + P_ee*nue_e + Pbar_ee*nuebar_e
                  + P_mumu*numu_misid_e + Pbar_mumu*numubar_misid_e + NC_e
Energies in GeV, events per bin."""
    write("NOvA", edges, table, header)


if __name__ == "__main__":
    print("DUNE");    build_dune()
    print("Hyper-K and T2K"); build_hyperk_and_t2k()
    print("NOvA");    build_nova()
