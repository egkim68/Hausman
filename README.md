Hausman Test Monte Carlo Simulation
Overview
This repository contains an R script for a Monte Carlo simulation studying the impact of missing data on the Hausman test's performance in panel data settings. It evaluates various panel architectures, model complexities, and missing data mechanisms.
Requirements

R ≥ 4.0.0
Packages: plm, dplyr, tidyr, purrr, ggplot2, scales, knitr, broom

Usage

Set your working directory in R.
Run: source("hausman_monte_carlo.R")
Results (tables, figures, summary) are saved to the ./results/ directory.

Simulation Design

Panel Architectures: Wide (N=400, T=4), Square (N=200, T=8), Long (N=100, T=16), each with 1,600 observations.
Model Complexities: Simple (1 variable), Standard (3), Complex (6), High-Dimensional (10).
Missing Data Mechanisms: Random (MCAR), Early Exit, Late Missing, Cyclical (MAR).
Dropout Rates: 10%, 20%, 30%, 40%.
Replications: 100 per scenario, across 144 feasible scenarios.

Outputs

Tables: Statistical comparisons, failure rates, specificity, ANOVA results, and pairwise comparisons.
Figures: Failure rate heatmap and specificity bar plot.
Files: Raw simulation data, failure analysis, excluded scenarios, and a summary report.

Notes

Runtime: ~15–30 minutes.
Set seed (1978) for reproducibility.
Results may vary slightly due to stochastic simulation.

License
MIT License
