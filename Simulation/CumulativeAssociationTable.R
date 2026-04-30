# Make simulation cumulative association table
# This needs to be run after final_simulation_0fill_revised.R

rm(list = ls())
gc()

library(readr)
library(dplyr)
library(tidyr)


for (min_gest_age in c(30, 37)) {
  # ---------------------------------------- #
  # ---------------------------------------- #
  # create combined final table that includes
  # all scenarios
  # ---------------------------------------- #
  # ---------------------------------------- #

  # find the files containing simulation results
  filenames <- list.files(
    path = paste0("~/ZeroFillDLM/Simulation/Output/mingest", min_gest_age, "/")
  )
  filenames <- filenames[grep("final_table_scenario", filenames)]

  # load and combine the tables
  combined_final_table <- NULL
  for (i in filenames) {
    tabletemp <- read_csv(paste0(
      "~/ZeroFillDLM/Simulation/Output/mingest",
      min_gest_age,
      "/",
      i
    ))
    tabletemp$file <- i
    combined_final_table <- bind_rows(combined_final_table, tabletemp)
  }

  # save results
  write_csv(
    combined_final_table,
    paste0(
      "~/ZeroFillDLM/Simulation/FinalTablesFigures",
      "/combined_final_table_",
      min_gest_age,
      ".csv"
    ),
  )
}
