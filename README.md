# Reproducing results from "Assessing the Impact of Covariate Distribution and Positivity Violation on Weighting-Based Indirect Comparisons: a Simulation Study" 

This repo allows to reproduce the results presented in the article: "Assessing the Impact of Covariate Distribution and Positivity Violation on Weighting-Based Indirect Comparisons: a Simulation Study, Serret-Larmande et al., 2025, DOI: xxx/xxxx". The source code is written in R version 4.5.1.

Steps to reproduce the results:
- With a working installation of renv, set up your R environment with 
- Source the `run_simulations.R` file, after adjusting the `options(mc.cores = N)` depending on your hardware
- The outut will be written in `results_simulations/yyyymmdd_hhmmss/`. 
