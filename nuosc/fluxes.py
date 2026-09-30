"""
Unoscillated neutrino fluxes at the far detectors
=================================================

All loaders return a `Flux` object with the *differential* flux

        dPhi/dE   in   neutrinos / cm^2 / GeV / POT

at the far detector, for the flavours 'numu', 'numubar', 'nue', 'nuebar'.
Each flavour keeps its own histogram (bin edges + values) because the
original files do not all use the same binning.

Data sources (see fluxes/README.md for details)
-----------------------------------------------
* DUNE : Optimized Engineered Nov 2017 beam, GLoBES flux files at the FD
         (1300 km), https://glaucus.crc.nd.edu/DUNEFluxes/
         units in file: nu / m^2 / GeV / POT, 0.25 GeV bins (bin centres given)
* T2K  : 2020 flux release, Super-Kamiokande (295 km, 2.5 deg off-axis),
         https://t2k-experiment.org/results/neutrino-beam-flux-prediction-2020/
         units in file: nu / cm^2 / 50 MeV / 1e21 POT, variable bins
* NOvA : 2017 flux release, NEAR detector only (FERMILAB-DATA-2017-01),
         https://publicdocs.fnal.gov/cgi-bin/ShowDocument?docid=8
         units in file: nu / m^2 / 1e6 POT per bin.
         The FD flux is *approximated* as ND x (L_ND / L_FD)^2, i.e. a point
         source. The real far/near ratio is not exactly 1/L^2 (the ND sees an
         extended source), so NOvA FD rates are only indicative (~10-20%).
"""

from pathlib import Path

from .paths import ROOT

import numpy as np

# Default location of the flux files: <project>/fluxes
FLUX_DIR = ROOT / "fluxes"

FLAVOURS = ("numu", "numubar", "nue", "nuebar")

# flavour name -> (index in nuosc: 0=e, 1=mu, 2=tau, is antineutrino)
FLAVOUR_INFO = {
    "nue": (0, False), "numu": (1, False), "nutau": (2, False),
    "nuebar": (0, True), "numubar": (1, True), "nutaubar": (2, True),
}

FLAVOUR_LATEX = {
    "nue": r"\nu_e", "numu": r"\nu_\mu", "nutau": r"\nu_\tau",
    "nuebar": r"\bar\nu_e", "numubar": r"\bar\nu_\mu", "nutaubar": r"\bar\nu_\tau",
}


class Flux:
    """
    Histogrammed differential flux dPhi/dE [nu / cm^2 / GeV / POT].

    Attributes
    ----------
    name, mode : e.g. 'DUNE', 'FHC'
    edges      : dict flavour -> bin edges [GeV]
    values     : dict flavour -> dPhi/dE in each bin
    """

    def __init__(self, name, mode, edges, values, note=""):
        self.name = name
        self.mode = mode
        self.edges = {f: np.asarray(e, float) for f, e in edges.items()}
        self.values = {f: np.asarray(v, float) for f, v in values.items()}
        self.note = note
        self.interpolate = True

    @property
    def flavours(self):
        return tuple(self.values)

    def __call__(self, flavour, E, interpolate=None):
        """
        dPhi/dE of `flavour` evaluated at energies E (0 outside the histogram).

        interpolate = True  : linear interpolation between bin centres (smooth,
                              default, see self.interpolate)
        interpolate = False : piecewise constant, i.e. the histogram itself
        """
        E = np.asarray(E, float)
        if flavour not in self.values:
            return np.zeros_like(E)
        edges, vals = self.edges[flavour], self.values[flavour]
        if interpolate is None:
            interpolate = self.interpolate
        if interpolate:
            centres = 0.5 * (edges[:-1] + edges[1:])
            out = np.interp(E, centres, vals)
            out[(E < edges[0]) | (E > edges[-1])] = 0.0
            return out
        idx = np.searchsorted(edges, E, side="right") - 1
        inside = (idx >= 0) & (idx < len(vals))
        out = np.zeros_like(E)
        out[inside] = vals[idx[inside]]
        return out

    def integral(self, flavour, E_min=0.0, E_max=np.inf):
        """Integrated flux [nu / cm^2 / POT] between E_min and E_max."""
        edges, vals = self.edges[flavour], self.values[flavour]
        lo = np.clip(edges[:-1], E_min, E_max)
        hi = np.clip(edges[1:], E_min, E_max)
        return np.sum(vals * (hi - lo))

    def __repr__(self):
        s = [f"Flux({self.name}, {self.mode})  [nu / cm^2 / GeV / POT]"]
        for f in self.flavours:
            e = self.edges[f]
            s.append(f"  {f:8s}: {len(e) - 1:4d} bins, {e[0]:.2f}-{e[-1]:.1f} GeV, "
                     f"integral = {self.integral(f):.3e} /cm^2/POT")
        if self.note:
            s.append("  note: " + self.note)
        return "\n".join(s)


# --------------------------------------------------------------------------- #
#  Loaders for the individual file formats                                    #
# --------------------------------------------------------------------------- #
def load_dune_globes(path, mode="FHC", name="DUNE", bin_width=0.25):
    """
    DUNE GLoBES flux file. Columns:
        E_centre  nue  numu  nutau  nuebar  numubar  nutaubar
    in nu / m^2 / GeV / POT at the far detector.
    """
    d = np.loadtxt(path)
    E = d[:, 0]
    edges = np.append(E - bin_width / 2, E[-1] + bin_width / 2)
    cols = {"nue": 1, "numu": 2, "nuebar": 4, "numubar": 5}
    values = {f: d[:, i] * 1e-4 for f, i in cols.items()}   # m^-2 -> cm^-2
    return Flux(name, mode, {f: edges for f in cols}, values)


def load_t2k_sk(path, mode="FHC", name="T2K"):
    """
    T2K 2020 text table for Super-K. Columns:
        bin#  E_lo  -  E_hi  numu  numubar  nue  nuebar
    in nu / cm^2 / 50 MeV / 1e21 POT.
    """
    d = np.genfromtxt(path, skip_header=3, usecols=(1, 3, 4, 5, 6, 7))
    edges = np.append(d[:, 0], d[-1, 1])
    cols = {"numu": 2, "numubar": 3, "nue": 4, "nuebar": 5}
    values = {f: d[:, i] / 0.05 / 1e21 for f, i in cols.items()}   # -> per GeV per POT
    return Flux(name, mode, {f: edges for f in cols}, values)


def load_nova_nd(directory, mode="FHC", name="NOvA", L_ND=1.0, L_FD=810.0):
    """
    NOvA 2017 near-detector tables, one file per flavour:
        {mode}_Flux_{flavour}_NOvA_ND_2017.txt   with columns  emin emax nus fe
    `nus` is in nu / m^2 / 1e6 POT per bin. Scaled to the FD with (L_ND/L_FD)^2.
    """
    directory = Path(directory)
    edges, values = {}, {}
    scale = 1e-4 / 1e6 * (L_ND / L_FD) ** 2          # m^-2 -> cm^-2, 1e6 POT -> POT, ND -> FD
    for f in FLAVOURS:
        d = np.loadtxt(directory / f"{mode}_Flux_{f}_NOvA_ND_2017.txt", comments="#")
        edges[f] = np.append(d[:, 0], d[-1, 1])
        values[f] = d[:, 2] / (d[:, 1] - d[:, 0]) * scale
    note = f"FD flux approximated as ND x ({L_ND} km / {L_FD} km)^2"
    return Flux(name, mode, edges, values, note=note)


LOADERS = {
    "dune_globes": load_dune_globes,
    "t2k_sk": load_t2k_sk,
    "nova_nd": load_nova_nd,
}


def load_flux(experiment, mode="FHC", flux_dir=None):
    """
    Load the far-detector flux of an experiment defined in nuosc.experiments,
    e.g. load_flux("DUNE", "RHC").
    """
    from .experiments import get_experiment

    cfg = get_experiment(experiment)
    if "flux" not in cfg:
        raise ValueError(f"No flux file configured for {cfg['name']}")
    fcfg = cfg["flux"]
    flux_dir = Path(flux_dir) if flux_dir is not None else FLUX_DIR
    path = flux_dir / fcfg[mode]
    loader = LOADERS[fcfg["format"]]
    kwargs = fcfg.get("kwargs", {})
    return loader(path, mode=mode, name=cfg["name"], **kwargs)
