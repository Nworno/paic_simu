library("DataExplorer")
library("data.table")
library("gtsummary")
library("tidyr")
library("dplyr")
library("ggplot2")
dir_akiki <- "//Psl-s-dfs/nas/SPEPS/BIOSPIM/UF Biostatistique/akiki"
dir_akiki2 <- "//Psl-s-dfs/nas/SPEPS/BIOSPIM/UF Biostatistique/akiki2/BIOSTATISTICIEN"

load(file.path(dir_akiki, "akiki.RData"))
load(file.path(dir_akiki2, "akiki_tb.RData"))

akiki1 <- as_tibble(akiki, .name_repair = "unique")
akiki2 <- as_tibble(akiki_tb, .name_repair = "unique")

labels_akiki <- sapply(akiki, \(x) attr(x, "label"))
names_akiki <- names(akiki)

configure_report(add_plot_correlation = FALSE, add_plot_missing = FALSE)
DataExplorer::create_report(akiki)
DataExplorer::create_report(akiki2)

# Patients count 

akiki |> 
  select(randomization_unique, inc_bras) |> 
  gtsummary::tbl_summary(by = inc_bras) 

akiki2 |>
  select(grp_trt) |> 
  tbl_summary()

#TODO: reprendre d'ici, comparer les résultats obtenus avec ceux de la publication
# RRT-free days

# Mortalité à J60
names_akiki[labels_akiki == "Vivant à J60"]
akiki |> 
  select(randomization_unique, all_of(names_akiki[labels_akiki == "Vivant à J60"])) |> 
  gtsummary::tbl_summary(by = randomization_unique)


akiki |> 
  select(randomization_unique, inc_bras) |> 
  gtsummary::tbl_summary(by = inc_bras)



# Vivant à J60 sans EER
# NULL
# Vivant à J60
# NULL
# Sevrage de l'EER
# Délai jusqu'au sevrage de l'EER (depuis J0)
# Délai jusqu'au sevrage de l'EER (depuis le début de l'EER)
# Mortality at day 28, numbers of mechanical ventilation–free and vasopressor-free days, length of stay in intensive care unit and in the hospital, and dependence on renal-replacement therapy at day 28 and 60 did not differ significantly between the two study groups 
