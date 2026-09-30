# Unoscillated far-detector fluxes

Public flux predictions (no oscillations) used by `nuosc.fluxes`.
Load them with `nuosc.load_flux("DUNE", "FHC")`; everything is converted to
**ν / cm² / GeV / POT at the far detector**.

FHC = forward horn current (neutrino-enhanced beam), RHC = reverse horn current
(antineutrino-enhanced beam).

## DUNE — `DUNE/`

| file | content |
|------|---------|
| `DUNE_FD_FHC_flux_globes.txt` | FD flux, neutrino mode |
| `DUNE_FD_RHC_flux_globes.txt` | FD flux, antineutrino mode |

* Source: <https://glaucus.crc.nd.edu/DUNEFluxes/OptimizedEngineeredNov2017/>
  (original names `histos_g4lbne_v3r5p4_QGSP_BERT_OptimizedEngineeredNov2017_{neutrino,antineutrino}_LBNEFD_globes_flux.txt`).
  Optimized 3-horn design, 2.2 m target, as used in the DUNE TDR.
* Format (GLoBES): `E_centre nue numu nutau nuebar numubar nutaubar`, 0.25 GeV bins.
* Units: ν / m² / GeV / POT at 1300 km.

## T2K — `T2K/`

`T2Kflux2020.tar.gz` is the original archive; `t2kflux_2020_public_release/` is its content
(see `README.pdf` there). We use the Super-K text tables:

| file | content |
|------|---------|
| `t2kflux_2020_plus250kA_nominal_sk.txt`  | SK flux, +250 kA (FHC), nominal beam |
| `t2kflux_2020_minus250kA_nominal_sk.txt` | SK flux, −250 kA (RHC), nominal beam |

* Source: <https://t2k-experiment.org/results/neutrino-beam-flux-prediction-2020/>
* Format: `bin#  E_lo - E_hi  numu  numubar  nue  nuebar`, 138 variable-width bins.
* Units: ν / cm² / 50 MeV / 10²¹ POT (per 50 MeV, whatever the bin width),
  at 295.3 km and 2.5° off-axis.
* `*_nd280.txt` are the near-detector fluxes, `*_runcond*` the fluxes for the actual
  run conditions, `t2kflux_2020_covariance.root` the flux covariance.

## NOvA — `NOvA/`

| file | content |
|------|---------|
| `{FHC,RHC}_Flux_{numu,numubar,nue,nuebar}_NOvA_ND_2017.txt` | ND flux per flavour and horn current |

* Source: <https://publicdocs.fnal.gov/cgi-bin/ShowDocument?docid=8>
  (FERMILAB-DATA-2017-01, DOI 10.15484/1959359).
* Format: `emin emax nus fe` (fe = fractional uncertainty); bins differ per flavour.
* Units: ν / m² / 10⁶ POT **per bin** at the **near detector** (~1 km).
* ⚠️ No public FD flux: `nuosc` approximates it as ND × (1 km / 810 km)², i.e. a point
  source. The true far/near ratio differs (the ND sees an extended source), so NOvA FD
  rates are only indicative.
