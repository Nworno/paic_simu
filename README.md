# Reproducing results from "Assessing the Impact of Covariate Distribution and Positivity Violation on Weighting-Based Indirect Comparisons: a Simulation Study" 

This repo allows to reproduce the results presented in the article: "Assessing the Impact of Covariate Distribution and Positivity Violation on Weighting-Based Indirect Comparisons: a Simulation Study, Serret-Larmande et al., 2025, DOI: xxx/xxxx". The source code is written in R version 4.5.1.

Steps to reproduce the results:

## 1. Clone the repository

```
git clone https://github.com/Nworno/paic_simu.git
```

## 2. Run the simulation

From an R console, restore the adequate package versions using the `renv` package.

```
install.packages("renv")
renv::restore()
```

Source the `run_simulations.R` file, after adjusting the `options(mc.cores = N)` on the very first line of that script, depending on your hardware.

```
Rscript run_simulations.R
```

## 3. Access the results

The Monte Carlo simulation output will be written in `results_simulations/yyyymmdd_hhmmss/`. Figures will be written in the subfolder `plots_publication`, and performance metrics numerical results in `processed_results/df_stats.rds`. 
