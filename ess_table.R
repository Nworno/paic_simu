library(data.table)
library(dplyr)
library(tidyr)
library(knitr)
library(kableExtra)
source("env_variables.R")

dir_results  <- file.path("results_simulations", DATE_EXPERIMENT)
dir_proc_ess <- file.path(dir_results, "processed_results")

joined_results <- readRDS(file.path(dir_proc_ess, "joined_results.rds"))

# ---------------------------------------------------------------------------
# Recode model names for display and aggregate ESS
# ---------------------------------------------------------------------------
model_display   <- c(ML = "PSW", MAIC_1 = "MAIC-1", MAIC_2 = "MAIC-2")
model_levels    <- c("PSW", "MAIC-1", "MAIC-2")
anchored_levels <- c("Anchored", "Unanchored")

ess_summary <- joined_results |>
  dplyr::filter(!is.na(Ess.x)) |>
  dplyr::mutate(Model = dplyr::recode(Model, !!!model_display)) |>
  dplyr::group_by(Population_parameters_num, Estimator_num, Model, Anchored) |>
  dplyr::summarize(
    mean_ess = mean(Ess.x, na.rm = TRUE),
    sd_ess   = sd(Ess.x,   na.rm = TRUE),
    n_iter   = sum(!is.na(Ess.x)),
    .groups  = "drop"
  ) |>
  dplyr::mutate(
    Model    = factor(Model,    levels = model_levels),
    Anchored = factor(Anchored, levels = anchored_levels),
    cell     = sprintf("%.1f (%.1f) [%d]", mean_ess, sd_ess, n_iter)
  )

# ---------------------------------------------------------------------------
# Pivot to wide: one column per (model x anchored) combination
# ---------------------------------------------------------------------------
ess_wide <- ess_summary |>
  dplyr::select(Population_parameters_num, Estimator_num, Model, Anchored, cell) |>
  tidyr::pivot_wider(names_from = c(Model, Anchored), values_from = cell,
                     names_sep = "___") |>
  dplyr::arrange(Population_parameters_num, Estimator_num) |>
  dplyr::rename(DGM = Population_parameters_num, Estimator = Estimator_num)

# Column order: PSW Anchored, PSW Unanchored, MAIC-1 Anchored, ...
data_cols <- as.vector(outer(model_levels, anchored_levels, paste, sep = "___"))

# ---------------------------------------------------------------------------
# Build LaTeX table with two-level header
# ---------------------------------------------------------------------------
header_top <- c(" " = 2, setNames(rep(2L, length(model_levels)), model_levels))

latex_table <- ess_wide |>
  dplyr::select(DGM, Estimator, all_of(data_cols)) |>
  knitr::kable(
    format    = "latex",
    booktabs  = TRUE,
    linesep   = "",
    col.names = c("DGM", "Estimator", rep(anchored_levels, length(model_levels))),
    caption   = paste("Effective sample size (Kish formula) per scenario and estimator.",
                      "Each cell reports mean (SD) [n iterations].")
  ) |>
  kableExtra::add_header_above(header = header_top) |>
  kableExtra::kable_styling(latex_options = c("hold_position", "scale_down"))

output_path <- file.path(dir_proc_ess, "ess_table.tex")
writeLines(latex_table, output_path)
cat("ESS table written to", output_path, "\n")
