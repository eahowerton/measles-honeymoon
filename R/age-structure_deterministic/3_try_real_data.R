library(dplyr)
library(reshape2)

folder = "data/WHO-data/"

births = read.csv(paste0(folder, "BirthRatePer1000downloadDec2024.csv")) %>%
  melt(c("Series.Name", "Series.Code", "Country.Name", "Country.Code"), variable.name = "year", value.name = "birth_rate") %>%
  mutate(year = as.integer(substr(year, 2, 5)))
vacc = read.csv(paste0(folder, "measlesVaccCoverFirstDose.csv"))
pop = read.csv(paste0(folder, "PopdownloadDec2024.csv")) %>%
  melt(c("X"), variable.name = "year", value.name = "pop") %>%
  mutate(year = as.integer(substr(year, 2, 5)))

case_study_country = "Samoa"



