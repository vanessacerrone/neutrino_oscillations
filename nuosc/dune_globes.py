"""
DUNE far-detector simulation from the official GLoBES configuration
====================================================================

This is a small, pure-Python re-implementation of what GLoBES does with the
DUNE TDR configuration published with

    DUNE Collaboration, "Experiment Simulation Configurations Approximating
    DUNE TDR", arXiv:2103.04797 (ancillary files)

stored in  data/DUNE_GLoBES_2103.04797/dune_globes/.

Compared with nuosc.rates (flux x linear cross section, efficiency 1, toy
Gaussian resolution), it uses the official DUNE ingredients:

  * flux      : flux/*.txt          nu / m^2 / GeV / POT at the far detector
  * xsec      : xsec/xsec_{cc,nc}.dat  GENIE 2.12 sigma/E on argon [1e-38 cm^2/GeV/nucleon]
                (columns: log10(E/GeV), nue, numu, nutau, nuebar, numubar, nutaubar)
  * smearing  : smr/*.txt           migration matrices true E -> reconstructed E
  * efficiency: eff/post_*.txt      selection efficiency in each reconstructed bin
  * channels & samples ("rules"), baseline, density, exposure: DUNE_GLoBES.glb

How an event rate is computed (the GLoBES recipe)
-------------------------------------------------
For each "channel" (flux flavour alpha -> detected flavour beta, CC or NC):

  1. true-energy sampling bins k (width dE_k, centre E_k):
         n_k = NORM * Phi_alpha(E_k) * sigma_beta(E_k) * P(alpha->beta; E_k) * dE_k
     with NORM = @norm * @power * @time * target_mass / L^2   (GLoBES convention)
  2. smear to reconstructed bins i:     m_i = sum_k S_ik n_k
  3. apply the selection efficiency:    N_i = eff_i * m_i

A "rule" (= an analysis sample, e.g. nue appearance in FHC) is the sum of
signal channels plus background channels, inside an energy window.
NC channels ("NOSC") are not oscillated (for 3 flavours the NC rate does not
change with oscillations). nu_tau cross sections are set to zero in this
configuration, as stated in the paper.
"""

import re
from pathlib import Path

from .paths import ROOT

import numpy as np

from .formulas import prob_matter_nufast, prob_vacuum

DUNE_GLOBES_DIR = ROOT / "data" / "DUNE_GLoBES_2103.04797" / "dune_globes"

_FLAV = {"e": 0, "m": 1, "t": 2}
_XSEC_COL = {("e", "+"): 1, ("m", "+"): 2, ("t", "+"): 3,
             ("e", "-"): 4, ("m", "-"): 5, ("t", "-"): 6}
_FLUX_COL = _XSEC_COL        # same column order in the GLoBES flux files


# --------------------------------------------------------------------------- #
#  small parsers for the GLoBES (AEDL) files                                  #
# --------------------------------------------------------------------------- #
def _strip_comments(text):
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return re.sub(r"//[^\n]*", "", text)


def _read(path):
    return _strip_comments(Path(path).read_text())


def _list(text, name):
    """Parse  $name = {a, b, c}  or  %name = {a, b, c}."""
    m = re.search(r"[\$%]" + name + r"\s*=\s*\{([^}]*)\}", text)
    return np.array([float(x) for x in m.group(1).replace("\n", " ").split(",") if x.strip()])


def _scalar(text, name):
    m = re.search(r"[\$@]?" + name + r"\s*=\s*([-+0-9.eE]+|[A-Z_]+)", text)
    return m.group(1)


def _parse_smearing(path, n_sampling):
    """GLoBES energy() block: one row {k_first, k_last, values...} per reconstructed bin."""
    rows = re.findall(r"\{([^}]*)\}", _read(path))
    S = np.zeros((len(rows), n_sampling))
    for i, row in enumerate(rows):
        v = [float(x) for x in row.split(",") if x.strip()]
        k0, k1 = int(v[0]), int(v[1])
        S[i, k0:k1 + 1] = v[2:2 + (k1 - k0 + 1)]
    return S


class DUNEGlobes:
    """
    DUNE far-detector event rates from the official GLoBES configuration.

    >>> dune = DUNEGlobes()
    >>> osc = NeutrinoOscillator("NO")
    >>> s = dune.sample(osc, "nue_app")     # FHC nu_e appearance sample
    >>> s["total"].sum(), s["signal"].sum()

    Samples ("rules"): 'nue_app' (FHC), 'nuebar_app' (RHC),
                       'numu_dis' (FHC), 'numubar_dis' (RHC).
    """

    def __init__(self, directory=None, years_fhc=None, years_rhc=None, mass_kt=None):
        self.dir = Path(directory) if directory is not None else DUNE_GLOBES_DIR
        glb = _read(self.dir / "DUNE_GLoBES.glb")
        defs = _read(self.dir / "definitions.inc")
        beam = _read(self.dir / "flux" / "Beam.inc")

        # ---- experiment set-up ---------------------------------------- #
        self.baseline = float(_list(glb, "lengthtab")[0])          # km
        self.rho = float(_list(glb, "densitytab")[0])               # g/cm^3
        self.years = {"FHC": float(_scalar(defs, "NUTIME")) if years_fhc is None else years_fhc,
                      "RHC": float(_scalar(defs, "NUBARTIME")) if years_rhc is None else years_rhc}
        self.mass_kt = float(_scalar(defs, "LAMASS")) if mass_kt is None else mass_kt
        self.power = float(re.search(r"@power\s*=\s*([0-9.eE+-]+)", beam).group(1))   # 1e20 POT / yr
        self.norm = float(re.search(r"@norm\s*=\s*([0-9.eE+-]+)", beam).group(1))

        # ---- binning -------------------------------------------------- #
        emin = float(_scalar(glb, r"\$emin"))
        self.sampling_edges = emin + np.concatenate([[0], np.cumsum(_list(glb, "sampling_stepsize"))])
        self.reco_edges = emin + np.concatenate([[0], np.cumsum(_list(glb, "binsize"))])
        self.E_true = 0.5 * (self.sampling_edges[1:] + self.sampling_edges[:-1])
        self.dE_true = np.diff(self.sampling_edges)
        self.E_reco = 0.5 * (self.reco_edges[1:] + self.reco_edges[:-1])

        # ---- flux and cross sections ---------------------------------- #
        self.flux = {}
        for mode, fname in re.findall(r"nuflux\(#flux_(\w+)\)<\s*@flux_file\s*=\s*\"([^\"]+)\"", beam):
            self.flux[mode] = np.loadtxt(self.dir / fname)
        self.xsec = {"CC": np.loadtxt(self.dir / "xsec" / "xsec_cc.dat"),
                     "NC": np.loadtxt(self.dir / "xsec" / "xsec_nc.dat")}

        # ---- efficiencies --------------------------------------------- #
        self.eff = {}
        for f in sorted((self.dir / "eff").glob("post_*.txt")):
            txt = _read(f)
            for name in re.findall(r"%(\w+)\s*=", txt):
                self.eff[name] = _list(txt, name)

        # ---- channels ------------------------------------------------- #
        pat = (r"channel\(#(\w+)\)<\s*@channel\s*=\s*#flux_(\w+)\s*:\s*([+-])\s*:\s*(\w)\s*:\s*"
               r"(NOSC_)?(\w)\s*:\s*#(\w+)\s*:\s*#(\w+)\s*@post_smearing_efficiencies\s*=\s*copy\(%(\w+)\)")
        self.channels = {}
        smear_cache = {}
        for name, mode, pol, fi, nosc, ff, xs, smr, eff in re.findall(pat, glb):
            if smr not in smear_cache:
                smear_cache[smr] = _parse_smearing(self.dir / "smr" / f"{smr}.txt", len(self.E_true))
            self.channels[name] = dict(mode=mode, polarity=pol, initial=fi, final=ff,
                                       oscillate=(nosc == ""), xsec=xs,
                                       smear=smear_cache[smr], eff=self.eff[eff])

        # ---- rules (analysis samples) --------------------------------- #
        self.rules = {}
        for name, body in re.findall(r"rule\(#(\w+)\)<(.*?)>", glb, flags=re.S):
            sig = re.search(r"@signal\s*=\s*([^\n]*)", body).group(1)
            bkg = re.search(r"@background\s*=\s*([^\n]*)", body).group(1)
            win = re.search(r"@energy_window\s*=\s*([0-9.]+)\s*:\s*([0-9.]+)", body)
            self.rules[name] = dict(signal=re.findall(r"#(\w+)", sig),
                                    background=re.findall(r"#(\w+)", bkg),
                                    window=(float(win.group(1)), float(win.group(2))))

    # ------------------------------------------------------------------ #
    def normalisation(self, mode):
        """GLoBES normalisation factor NORM = norm * power * time * mass / L^2."""
        return self.norm * self.power * self.years[mode] * self.mass_kt / self.baseline ** 2

    def pot(self, mode):
        return self.power * 1e20 * self.years[mode]

    def flux_at(self, mode, flavour, polarity, E):
        """Unoscillated flux [nu / m^2 / GeV / POT] (linear interpolation, as GLoBES)."""
        d = self.flux[mode]
        return np.interp(E, d[:, 0], d[:, _FLUX_COL[(flavour, polarity)]], left=0, right=0)

    def xsec_at(self, kind, flavour, polarity, E):
        """Cross section [1e-38 cm^2 per nucleon] = (sigma/E)(E) * E, interpolated in log10(E)."""
        d = self.xsec[kind]
        s_over_E = np.interp(np.log10(E), d[:, 0], d[:, _XSEC_COL[(flavour, polarity)]], left=0, right=0)
        return s_over_E * E

    def probability_matrix(self, osc, polarity, matter=True):
        """Full 3x3 probability matrix at the true sampling-bin centres."""
        anti = polarity == "-"
        if matter:
            return prob_matter_nufast(self.E_true, self.baseline, **osc.params, rho=self.rho,
                                      Ye=osc.Ye, antineutrino=anti)
        return prob_vacuum(self.E_true, self.baseline, **osc.params, antineutrino=anti)

    def _probabilities(self, osc, matter):
        if osc is None:
            return None
        return {pol: self.probability_matrix(osc, pol, matter) for pol in "+-"}

    def channel_true(self, osc, name, matter=True, _P=None):
        """Events per true-energy sampling bin (before smearing and efficiency)."""
        ch = self.channels[name]
        E = self.E_true
        phi = self.flux_at(ch["mode"], ch["initial"], ch["polarity"], E)
        xs = self.xsec_at(ch["xsec"], ch["final"], ch["polarity"], E)
        if not ch["oscillate"]:
            P = 1.0                                   # NC: unaffected by 3-flavour oscillations
        elif osc is None:
            P = 1.0 if ch["initial"] == ch["final"] else 0.0
        else:
            Pm = _P if _P is not None else self._probabilities(osc, matter)
            P = Pm[ch["polarity"]][_FLAV[ch["initial"]], _FLAV[ch["final"]]]
        return self.normalisation(ch["mode"]) * phi * xs * P * self.dE_true

    def channel(self, osc, name, matter=True, _P=None):
        """Selected events per reconstructed-energy bin for one channel."""
        ch = self.channels[name]
        return ch["eff"] * (ch["smear"] @ self.channel_true(osc, name, matter, _P))

    def sample(self, osc, rule, matter=True, window=True, _P=None):
        """
        All channels of an analysis sample, in reconstructed-energy bins.

        Returns a dict with every channel, plus 'signal', 'background', 'total',
        'edges' and 'mask' (reco bins inside the rule's energy window).
        osc = None -> no oscillations.
        """
        if _P is None:
            _P = self._probabilities(osc, matter)
        r = self.rules[rule]
        out = {ch: self.channel(osc, ch, matter, _P) for ch in r["signal"] + r["background"]}
        out["signal"] = sum(out[ch] for ch in r["signal"])
        out["background"] = sum(out[ch] for ch in r["background"])
        out["total"] = out["signal"] + out["background"]
        lo, hi = r["window"]
        out["mask"] = (self.reco_edges[:-1] >= lo - 1e-9) & (self.reco_edges[1:] <= hi + 1e-9)
        if window:
            for k in list(out):
                if k != "mask":
                    out[k] = np.where(out["mask"], out[k], 0.0)
        out["edges"] = self.reco_edges
        return out

    def all_samples(self, osc, matter=True, window=True):
        """dict rule -> sample(...) for the four DUNE samples (probabilities computed once)."""
        P = self._probabilities(osc, matter)
        return {r: self.sample(osc, r, matter, window, _P=P) for r in self.rules}

    def __repr__(self):
        return (f"DUNEGlobes(L = {self.baseline} km, rho = {self.rho} g/cm^3, "
                f"{self.mass_kt:g} kt, {self.years['FHC']:g} yr FHC + {self.years['RHC']:g} yr RHC "
                f"at {self.power:g}e20 POT/yr)\n"
                f"  {len(self.E_true)} true sampling bins, {len(self.E_reco)} reco bins, "
                f"{len(self.channels)} channels, samples: {list(self.rules)}")
