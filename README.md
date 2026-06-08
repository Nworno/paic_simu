# Reproducing results from "Assessing the Impact of Covariate Distribution and Positivity Violation on Weighting-Based Indirect Comparisons: a Simulation Study" 

Theses files allows to reproduce the results presented in the article: "Assessing the Impact of Covariate Distribution and Positivity Violation on Weighting-Based Indirect Comparisons: a Simulation Study, Serret-Larmande et al.". The source code is written in R version 4.5.1.

Steps to reproduce the results:

# Run the simulation

From an R console, restore the adequate package versions using the `renv` package.

```
install.packages("renv")
renv::activate()
renv::restore()
```

Source the `run_simulations.R` file, after adjusting the number of iterations and the number of bootstrap replications.

```
Rscript run_simulations.R
```

## Access the results

The Monte Carlo simulation output will be written in `results_simulations/yyyymmdd_hhmmss/`. Figures will be written in the subfolder `plots_publication`, and performance metrics numerical results in `processed_results/df_stats.rds`. 
