# Hausman Test Monte Carlo Reproducibility

Code accompanying **"Wider Panels, Higher Hausman Test Specificity: A Monte Carlo Study Under Missing Data"** (Eungi Kim).

This repository contains two independent implementations of the same 192,000-run Monte Carlo design (192 scenarios × 1,000 replications: 3 panel architectures × 4 model complexity levels × 4 missing-data mechanisms × 4 missingness rates). They share no code and were built directly from the manuscript's equations, so results from one can be checked against the other.

## Files

- `hausman_full_reproducibility.R` — primary implementation (R, `plm`). Regenerates the raw simulation results and every table (1–5) reported in the manuscript.
- `hausman_full_reproducibility.py` — independent implementation (Python, `pandas` + `linearmodels`). Regenerates the raw simulation results used for the cross-validation reported in Section 5.5.

## Requirements

**R (4.3.3):** `plm`, `dplyr`, `tidyr`, `purrr`, `car`, `ggplot2`

**Python:** `numpy`, `pandas`, `linearmodels`, `scipy`

## Running

```r
# R
source("hausman_full_reproducibility.R")
```

```bash
# Python
python hausman_full_reproducibility.py
```

Both scripts write to a local `results/` folder and checkpoint after every scenario, so an interrupted run can be resumed by simply rerunning the script. No external data are used — all panels are simulated. Seed is fixed (`1978`) for exact reproducibility; minor numerical differences may occur across software versions.

## Output

- Raw replication-level results (one row per simulation run)
- Table CSVs corresponding to the manuscript's Tables 1–5

## Citation

If you use this code, please cite the manuscript above.
