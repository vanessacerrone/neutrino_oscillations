# Data Release 2024

## Introduction

This data release corresponds to the Bayesian 2024 analysis of NOvA $\nu_e$ appearance and $\nu_{\mu}$ disappearance data, corresponding to analysis described in [https://doi.org/10.1103/x53y-2b86](https://doi.org/10.1103/x53y-2b86).

File `NOvA_2024_data_histograms.root` contains data histograms for all the NOvA data samples.

Exposure:

* neutrino-enhanced beam: 26.61E20 protons on target
* antineutrino-enhanced beam: 12.5E20 protons on target

External constraints:

* ss2th12=0.851, dm21=7.53e-5 are fixed at 2019 PDG values, with negligible effect on NOvA predictions.
* ss2th13 & dm32:
  * RCDB1D: 1D constraint from Daya Bay for ss2th13, 0.0851+/-0.0024
  * RCDB2D: Correlated 2D constraint from Daya Bay on ss2th13 & dm32, available in their official data release: [https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.130.161802#supplemental](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.130.161802#supplemental)

An additional file containing the predictions (`NOvA_2024_prediction_with_systs_histograms.root`) for all channels’ signal and background components computed at the NOvA best-fit oscillation parameters and systematic pull terms. Numu samples include the no-oscillation case as well.

Please note that any fits performed with these histograms are not expected to exactly reproduce the official NOvA results as parameterizations of the numerous systematic uncertainties considered in the official fits are not included in this release.

## Credible intervals

Root files that begin with 'contours_' contain 2D credible interval contours in `TGraph` objects, with `axis_` objects corresponding to the 2D axis for the credible intervals.

Key for decoding the file names and objects in these files:

* **Filenames** are formatted as `contours__{NOvAOnly,RCDB1D,RCDB2D}_{conditional,nonconditional}.root`
  * `NOvAOnly`: No external oscillation constraints applied.
  * `RCDB1D`: 1D Reactor Constraint from DayaBay 2024 results applied to $\sin^2 2\theta_{13}$.
  * `RCDB2D`: 2D Reactor Constraint from DayaBay 2024 results applied to $\sin^2 2\theta_{13}$ and $\Delta m^2_{32}$.
  * `conditional`: Marginalized over each mass ordering separately.
  * `nonconditional`: Marginalized over both mass orderings simultaneously.
* **Objects in the files** are named with a long common prefix, and separate suffices in a format: `<variable><variable>_{NO|IO|BO}_[cred_int_{068|06827|090|095|09545|099|09973}]`:
  * `<variable><variable>` are the variable name:
    - `dcp`: $\delta_{CP}$
    - `sindcp`: $\sin \delta_{CP}$
    - `ss2th13`: $\sin^2 2\theta_{13}$
    - `dm32`: $\Delta m^2_{32}$
    - `ssth23`: $\sin^2 \theta_{23}$
    - `J`: Jarlskog invariant
  * `NO`: Normal mass ordering
  * `IO`: Inverted mass ordering
  * `BO`: Both mass orderings
  * The number following `cred_int`, if present, indicates the fraction of the posterior contained by the credible region: 68% 68.27%, 90%, 95%, 95.45%, 99%, 99.73%
  * TGraphs corresponding to contours also contain a final numerical suffix differentiating disconnected components at the same contour level (All of the TGraphs that differ only by this final suffix should be plotted together to obtain the final contour).
