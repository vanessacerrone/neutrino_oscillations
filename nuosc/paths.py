"""
Where the data folders (fluxes/, data/, far_detector_spectra/) live
===================================================================

The data are not part of the Python package: they stay in the project folder.
The project folder is found, in this order, as

  1. the environment variable NUOSC_ROOT, if set;
  2. the folder that contains the nuosc/ package (editable install, `pip install -e .`,
     or running from the project folder);
  3. the current working directory or one of its parents.

Only the probability code (NeutrinoOscillator, prob_vacuum, prob_matter_nufast, ...)
works without the data folders.
"""

import os
from pathlib import Path

_MARKERS = ("far_detector_spectra", "fluxes", "data")


def _looks_like_root(p):
    return any((p / m).is_dir() for m in _MARKERS)


def project_root():
    env = os.environ.get("NUOSC_ROOT")
    if env:
        return Path(env).expanduser().resolve()
    here = Path(__file__).resolve().parent.parent
    if _looks_like_root(here):
        return here
    cwd = Path.cwd().resolve()
    for p in (cwd, *cwd.parents):
        if _looks_like_root(p):
            return p
    return here            # fallback: loaders will raise a clear FileNotFoundError


ROOT = project_root()
