# Spectral‑domain test for spatial variation in periodically correlated hydrological time series

This repository contains the R implementation of the method described in the paper:

> "Spatial variation analysis of periodically correlated hydrological time series: A new insight"

## Contents

- `PCspatial_analysis.R` – Main R script with all core functions and simulation examples
- `Data_spatial_spectral.csv` – Processed monthly precipitation (9 stations) and runoff (4 stations) data from Golestan Province, Iran, in CSV format.


## Requirements

- R (version 4.3.0 or higher)
- Required packages: `MASS`, `dplyr`, `tidyr`, `ggplot2` (optional for plotting)

Install the required packages in R:
```r
install.packages(c("MASS", "dplyr", "tidyr", "ggplot2"))