# Finite-Density Compactness Threshold Framework

This repository contains the reproducibility code, validation outputs, and figure-generation pipeline for the paper on finite-density compactness thresholds, EOS--TOV reconstruction, and tidal-deformability estimates.

## Contents

- `scripts/01_download_compose_eos.py`: downloads APR, SLY4, and DD2 EOS tables from CompOSE.
- `scripts/02_threshold_scan.py`: computes the finite-density compactness-threshold mass scan.
- `scripts/03_full_tov_reconstruction.py`: performs full thermo--TOV reconstruction from CompOSE thermodynamic tables and validates against official `eos.mr` sequences.
- `scripts/04_generate_figures.py`: generates the mass--radius, compactness, thermo--TOV validation, and approximate tidal-deformability figures.
- `outputs/figures/`: final paper figures.
- `outputs/tables/`: CSV tables from threshold scans and EOS comparisons.
- `outputs/validation/`: thermo--TOV validation summaries and reconstructed mass--radius sequences.

## EOS Models

APR, SLY4, and DD2 from the CompOSE database.

## Main Validation Result

| EOS | Mass error | Radius error |
|---|---:|---:|
| APR | 0.18% | 0.04% |
| SLY4 | 0.13% | 0.19% |
| DD2 | 0.32% | 1.48% |

## Run Order

```bash
pip install -r requirements.txt
python scripts/01_download_compose_eos.py
python scripts/02_threshold_scan.py
python scripts/03_full_tov_reconstruction.py
python scripts/04_generate_figures.py
```

## Notes

The thermo--TOV reconstruction is performed from the CompOSE `eos.nb` and `eos.thermo` thermodynamic tables and validated against the official CompOSE `eos.mr` mass--radius sequences. The tidal-deformability estimates are approximate compactness-based diagnostics and are not a substitute for full relativistic Love-number integration.

## Citation

If using this repository, cite the corresponding Zenodo release.
