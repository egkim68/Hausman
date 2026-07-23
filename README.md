# Hausman Test Under Missing Data: Monte Carlo Simulation

R code for the Monte Carlo study in "The Impact of Missing Data on the Hausman Test: A Monte Carlo Study" (Eungi Kim).

## Requirements

R (>= 4.0.0) with packages: `plm`, `dplyr`, `tidyr`, `purrr`, `ggplot2`, `scales`, `knitr`, `broom`

## Usage

```r
source("hausman_monte_carlo.R")
```

Running the script regenerates all 14,400 simulation runs and produces every table and figure reported in the manuscript, saved to a `results/` folder created automatically at runtime. `set.seed(1978)` is fixed at the top of the script for determinism within a given R/package version.

## Note on exact reproducibility

Because results depend on R's random-number generation, minor numerical differences (e.g., a handful of individual run outcomes near a p = 0.05 threshold) are possible across different R or package versions. The qualitative findings (which missing-data mechanisms fail vs. succeed, the panel architecture/complexity ordering) are robust to this.

## Contact

Eungi Kim, Department of Library and Information Science, Keimyung University, Daegu, Republic of Korea
egkim@kmu.ac.kr

## License

MIT
