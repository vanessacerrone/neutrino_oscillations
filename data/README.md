# Public far-detector data

## NOvA 2024 — `NOvA_2024_FD_release/`

* Source: <https://novaexperiment.fnal.gov/data-releases/> → "Oscillation results using 2024
  Bayesian approach", <https://publicdocs.fnal.gov/cgi-bin/ShowDocument?docid=594>
  (FERMILAB-DATA-2025-05, DOI 10.5281/zenodo.17822358).
* Exposure: 26.61×10²⁰ POT neutrino beam (FHC) + 12.5×10²⁰ POT antineutrino beam (RHC).
* `NOvA_2024_data_release.zip`: original archive; `data_release_2024/`: its content
  (ROOT files + NOvA's `README.md`).
* `csv/`: plain-text versions made from the ROOT files (no ROOT/uproot needed), read by `nuosc.data`:

| file | content |
|------|---------|
| `NOvA_2024_FD_numu_FHC.csv`, `NOvA_2024_FD_numu_RHC.csv` | νμ CC sample vs reconstructed energy (all quartiles summed): data, NOvA best-fit prediction (total, signal, beam bkg), no-oscillation prediction, cosmics |
| `NOvA_2024_FD_nue_totals.csv` | νe samples, total events: data and prediction by component (the νe histograms use PID × energy "analysis bins") |

The NOvA predictions are at NOvA's best-fit oscillation parameters and systematic pulls; fits to
these histograms will not exactly reproduce the official NOvA results (see NOvA's README).

## DUNE TDR GLoBES configuration — `DUNE_GLoBES_2103.04797/`

* Source: ancillary files of DUNE Collaboration, *Experiment Simulation Configurations Approximating
  DUNE TDR*, [arXiv:2103.04797](https://arxiv.org/abs/2103.04797).
* `dune_globes/`: the GLoBES configuration (`DUNE_GLoBES.glb`, `definitions.inc`, flux, GENIE cross
  sections, smearing matrices, efficiencies), read directly by `nuosc.dune_globes.DUNEGlobes`.
* `dune_flux/`: ND and FD fluxes (ROOT + GLoBES text).
* `*.png`, `tdr_configs_arxiv.tex`: figures and source of the paper, used in
  `dune_tdr_simulation.ipynb` to validate our implementation.
* Default exposure: 6.5 + 6.5 years at 1.1e21 POT/yr (1.2 MW), 40 kt = 624 kt-MW-yr ("10 years staged").

## T2K Run 1-9 — `T2K_OA2019_release/`

* Source: <https://t2k-experiment.org/results/t2kdata-oa-2019/> (data release of
  [arXiv:2101.03779](https://arxiv.org/abs/2101.03779), PRD 103, 112008), 14.94e20 POT FHC + 16.35e20 POT RHC.
* `Analysis_{A,B,C}.root`, `README`: official release (confidence regions, Δχ², posteriors).
* T2K does not release far-detector spectra. `csv/` contains, as plain text:

| file | content |
|------|---------|
| `T2K_OA2019_SK_event_rates.csv` | Table XI of the paper: predicted events per SK sample for δCP = −π/2, 0, π/2, π, and observed events |
| `T2K_OA2019_A_dchi2_dcp_wRC_{NO,IO}.csv` | official Δχ² vs δCP (analysis A, with reactor constraint), from `Analysis_A.root` |

* `T2K_2101.03779.pdf` (the paper, for reference) is ignored by git.

## Hyper-K Design Report — `HyperK_DR_1805.04163/`

* `HyperK_DesignReport_1805.04163.pdf`: Hyper-Kamiokande Design Report, arXiv:1805.04163 (ignored by git).
  Tables XXXVI–XXXVII (expected νe/νμ candidate events) are used to calibrate
  `far_detector_spectra/HyperK_FD_spectra.csv`; the numbers are copied in `tools/build_fd_spectra.py`.
