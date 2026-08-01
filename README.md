# Spatial variation analysis of periodically correlated hydrological time series


**Spatial variation analysis of periodically correlated hydrological time series: A new insight**

This repository provides the processed 15-year precipitation and runoff datasets, the original source workbook, R scripts for every analysis reported in the manuscript.



## Data

Both analysis datasets contain 180 monthly observations (15 years):

- `data/precipitation_monthly_180.csv`: nine meteorological stations.
- `data/runoff_monthly_180.csv`: four hydrometric stations.


## Software

Recommended: R 4.3 or later.

Install the package required:

```bash
Rscript R/00_install_packages.R
```

The spectral test, simulations, figures, and data checks otherwise use base R.

## Execution

Run commands from the repository root.


### Real-Data Reproduction

```bash
Rscript run_all.R real
```

This runs the precipitation and runoff spectral tests, Holm-adjusted pairwise tests, regional precipitation comparisons, descriptive summaries, figures, Friedman-Conover-Holm analyses, agreement calculations, and validation files using the full real-data permutation setting.

### Complete manuscript simulation design

```bash
Rscript run_all.R full
```

This runs 1,000 Monte Carlo replications per simulation cell and 2,000 within-frequency-block permutations per test.

## Analysis-to-script map

| Manuscript component | Script | Main output |
|---|---|---|
| Tables 2-7, Scenarios 1-3 | `R/02_simulation_scenarios_1_3.R` | `simulation_power_scenarios_1_3.csv` |
| Tables 8-9, Scenario 4 | `R/03_simulation_scenario_4.R` | `simulation_power_scenario_4.csv` |
| Scenario 4 scaled t3 innovations | `R/03_simulation_scenario_4.R` | `scenario_4_nongaussian_sensitivity.csv` |
| Empirical Type I error | `R/04_type_I_error.R` | `empirical_type_I_error.csv` |
| Daniell bandwidth sensitivity | `R/05_bandwidth_sensitivity.R` | `bandwidth_sensitivity_scenario_3.csv` |
| Real-data global and pairwise tests | `R/06_real_data_spectral_analysis.R` | spectral global/pairwise CSV files |
| West/Central/East precipitation contrasts | `R/06_real_data_spectral_analysis.R` | zone global/pairwise CSV files |
| Figures 3-4 | `R/06_real_data_spectral_analysis.R` | PNG files in `figures/generated/` |
| Friedman-Conover-Holm comparison | `R/07_conventional_friedman_conover.R` | conventional pairwise/summary CSV files |
| Input-data integrity | `R/08_data_checks.R` | `input_data_integrity_checks.csv` |


## Method implementation

The core script follows the manuscript formulation:

1. Compute the `2T`-dimensional real/imaginary DFT vector at each selected evaluation frequency.
2. Construct the block-specific spatial and residual periodogram-type quadratic forms corresponding to equations (6)-(7).
3. Apply a symmetric Daniell moving average to these quadratic-form blocks. Boundary windows are truncated and the remaining equal weights are renormalized.
4. Calculate Wilks' Lambda from the weighted matrices.
5. Generate the null distribution by independently permuting station labels within each frequency block.
6. Apply Holm adjustment to the family of pairwise permutation p-values.

The frequency grid is

```text
lambda_j = 2*pi*(j-1)/(k*T),  j = 1,...,k.
```

The adaptive Daniell half-width is

```text
m = min(floor(sqrt(n)), floor((k-1)/2)).
```

## Random seeds and output folders

All seeds and analysis settings are centralized in `R/00_config.R`. Newly generated files are written to:

- `results/generated/`
- `figures/generated/`

