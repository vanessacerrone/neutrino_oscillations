"""
Plot style (optional)
=====================

    from nuosc.plotting import set_style
    c = set_style(True)     # use the nuosc style; returns the list of colours
    c = set_style(False)    # change NOTHING: keep your own matplotlib settings;
                            # returns the colours of your current colour cycle

With enabled=True the style is chosen in this order:
  1. the matplotlib style named `style` (default "latex_style"), if you have it installed;
  2. SciencePlots ("science" + "no-latex"; installed with the package);
  3. the matplotlib default style;
and the IBM colour-blind-safe colour cycle is used.
"""

import matplotlib.pyplot as plt

IBM_COLOR_CYCLE = [
    '#648FFF', '#DC267F', '#FFB000', '#785EF0', 'forestgreen', '#FE6100',
    '#40B0A6', '#006CD1', '#B00068', '#D35FB7', '#994F00'
]


def current_colors():
    """Colours of the current matplotlib colour cycle."""
    return list(plt.rcParams["axes.prop_cycle"].by_key().get("color", IBM_COLOR_CYCLE))


def set_style(enabled=True, style="latex_style", colors=IBM_COLOR_CYCLE, verbose=True):
    """
    enabled : True  -> set the style described above;
              False -> do not touch matplotlib at all.
    style   : name (or list of names) of the matplotlib style to try first.
    colors  : colour cycle to use when enabled (None = keep the style's own cycle).
    Returns the list of colours of the (new or current) colour cycle.
    """
    if not enabled:
        if verbose:
            print("plot style: unchanged (your own matplotlib settings)")
        return current_colors()

    used = None
    try:
        plt.style.use(style)
        used = str(style)
    except (OSError, ValueError):
        try:
            import scienceplots  # noqa: F401  (registers the "science" styles)
            plt.style.use(["science", "no-latex"])
            plt.rcParams.update({"figure.dpi": 110})
            used = "scienceplots (science, no-latex)"
        except (ImportError, OSError, ValueError):
            plt.style.use("default")
            used = "matplotlib default"
    if colors is not None:
        plt.rcParams["axes.prop_cycle"] = plt.cycler("color", list(colors))
    if verbose:
        print(f"plot style: {used}")
    return current_colors()
