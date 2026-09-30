# Far-detector spectra for T2K, NOvA, DUNE and Hyper-K

One CSV per experiment, `<EXP>_FD_spectra.csv`, ready to be multiplied by oscillation
probabilities. Each file has a comment header (source, exposure, baseline, density, recipe)
followed by one row per energy bin:

| column | meaning (events per bin **if P = 1**) | multiply by |
|--------|----------------------------------------|-------------|
| `E_lo, E_hi, E_center` | bin edges and centre [GeV] (true neutrino energy) | |
| `<mode>_numu_mu` | ν_μ flux → ν_μ CC, μ-like sample | P(ν_μ→ν_μ) |
| `<mode>_numubar_mu` | ν̄_μ flux → ν̄_μ CC, μ-like sample | P(ν̄_μ→ν̄_μ) |
| `<mode>_numu_e` | ν_μ flux → ν_e CC, e-like sample (appearance signal) | P(ν_μ→ν_e) |
| `<mode>_numubar_e` | ν̄_μ flux → ν̄_e CC, e-like sample | P(ν̄_μ→ν̄_e) |
| `<mode>_nue_e`, `<mode>_nuebar_e` | intrinsic beam ν_e, ν̄_e → e-like sample | P(ν_e→ν_e), P(ν̄_e→ν̄_e) |
| `<mode>_numu_misid_e`, `<mode>_numubar_misid_e` | ν_μ CC mis-identified as e-like (DUNE only, 0 elsewhere) | P(ν_μ→ν_μ), P(ν̄_μ→ν̄_μ) |
| `<mode>_NC_mu`, `<mode>_NC_e` | neutral currents (+ cosmics for NOvA) | 1 |

`<mode>` = `FHC` (neutrino beam) or `RHC` (antineutrino beam). The prediction is

```
mu_like = P_mumu*numu_mu + Pbar_mumu*numubar_mu + NC_mu
e_like  = P_mue*numu_e + Pbar_mue*numubar_e + P_ee*nue_e + Pbar_ee*nuebar_e
          + P_mumu*numu_misid_e + Pbar_mumu*numubar_misid_e + NC_e
```

with probabilities evaluated at `E_center` (use the baseline and density in the header),
or simply `nuosc.predict("DUNE", NeutrinoOscillator("NO"))`. See `fd_spectra_quickstart.ipynb`.

## How each file was made (`python tools/build_fd_spectra.py`)

| experiment | ingredients | exposure | validation |
|------------|-------------|----------|------------|
| **DUNE** | official TDR GLoBES configuration (arXiv:2103.04797): flux, GENIE xsec on Ar, smearing matrices, efficiencies; selected events projected back on true energy (99 bins, 0–110 GeV) | 6.5 + 6.5 yr, 1.2 MW, 40 kt (624 kt-MW-yr) | totals identical to the full GLoBES calculation (`nuosc.dune_globes`) |
| **Hyper-K** | T2K 2020 SK flux (same beam, 2.5° off-axis, 295 km), GENIE σ/E (smoothed), one constant efficiency per component calibrated to Tables XXXVI–XXXVII of the Hyper-K Design Report (arXiv:1805.04163); e-like only for E < 1.25 GeV | 2.7e22 POT, ν:ν̄ = 1:3, 187 kt | reproduces the Design Report totals at its parameters |
| **T2K** | as Hyper-K, scaled to 22.5 kt and Run 1–9 POT, plus one normalisation per sample from Table XI of arXiv:2101.03779 | 14.94e20 + 16.35e20 POT | δCP dependence of Table XI reproduced to ≲ 4% |
| **NOvA** | μ-like: NOvA's own no-oscillation FD prediction (2024 release, reco energy); e-like: same FD ν_μ spectrum shape, normalised to Table 2 of arXiv:2509.04361 at the NOvA best fit | 26.61e20 + 12.5e20 POT | μ-like total 407.8 vs NOvA 408.6 at the NOvA best fit |

**Limitations:** energies are *true* energies (DUNE, T2K, Hyper-K) or NOvA's reconstructed energy
used as a proxy; T2K/Hyper-K efficiencies are energy independent and the e-like cut is applied in
true energy; NC shapes for T2K/Hyper-K/NOvA e-like are approximate; ν_τ CC events are not included.
For DUNE in reconstructed energy use `nuosc.dune_globes.DUNEGlobes`.
