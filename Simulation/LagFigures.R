# Make simulation cumulative association table
# This needs to be run after final_simulation_0fill_revised.R

rm(list = ls())
gc()

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)


# true exposure effect
truth <- rep(0, 42)
truth[15:30] <- c(1 / (1 + exp(seq(3, -4))), 1 / (1 + exp(seq(-4, 3))))
truth <- -truth / sum(truth) # rescale to have cumulative effect of -1
# plot(truth)

# this is used for plotting later
truth_df <- data.frame(week = 1:42, value = truth)


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
  filenames <- filenames[grep("sim_results_revised_scenario", filenames)]

  # load and combine the tables
  combined_results_all <- NULL
  for (i in filenames) {
    load(
      file = paste0(
        "~/ZeroFillDLM/Simulation/Output/mingest",
        min_gest_age,
        "/",
        i
      )
    )
    results_all$file <- i
    results_all$GAcontrol <- ifelse(
      grepl("GAcontrolTRUE", i),
      "adjusted for\ngestational age",
      "not adjusted for\ngestational age"
    )
    results_all$scenario <- ifelse(
      substring(i, 30, 30) == "B",
      "A) gestational age as predictor",
      ifelse(
        substring(i, 30, 30) == "C",
        "B) gestational age as mediator",
        ifelse(
          substring(i, 30, 30) == "E",
          "C) gestational age as predictor w/ confounding",
          "error"
        )
      )
    )

    combined_results_all <- bind_rows(combined_results_all, results_all)
  }

  #-------------------------------------------
  # plots by lag week
  #-------------------------------------------
  lag_plt_data_all <-
    combined_results_all |>
    select(-starts_with("cover")) |>
    filter(model == "ns (best AIC)") |>
    pivot_longer(
      cols = paste0("est", 1:42),
      names_to = "week",
      values_to = "est"
    ) |>
    mutate(week = as.numeric(substr(week, 4, 5))) |>
    select(-c("df", "AIC", "min_aic")) |>
    group_by(model, data, GAcontrol, scenario, week) |>
    summarize(
      mean = mean(est, na.rm = TRUE),
      min = min(est, na.rm = TRUE),
      max = max(est, na.rm = TRUE),
      q95 = quantile(est, 0.95, na.rm = TRUE),
      q05 = quantile(est, 0.05, na.rm = TRUE)
    )
  # this will produce warnings because of the 37 week group with all NAs after week 37

  # remove extra weeks from 37 week truncation so plot doesn't
  # look like the estimates goes to 0 in week 38
  lag_plt_data_all2 <- lag_plt_data_all |>
    filter(
      data != paste0("truncated at week ", min_gest_age) |
        week <= min_gest_age
    )

  # make plot with ns only for publication main text
  lag_plt_all_ns <-
    ggplot(
      data = lag_plt_data_all2 |>
        filter(model %in% c("ns (best AIC)"))
    ) +
    geom_ribbon(
      aes(x = week, y = mean, ymin = min, ymax = max),
      fill = "grey80"
    ) +
    geom_ribbon(
      aes(x = week, y = mean, ymin = q05, ymax = q95),
      fill = "grey50"
    ) +
    geom_line(aes(x = week, y = mean)) +
    # facet_grid(scenario+GAcontrol ~ data) +
    # facet_grid2(
    #   rows = vars(scenario, GAcontrol),
    #   cols = vars(data),
    #   strip = strip_nested()
    # ) +
    theme_bw(base_size = 9) +
    theme(strip.background = element_blank()) +
    facet_nested(scenario + GAcontrol ~ data, nest_line = TRUE) +
    geom_line(data = truth_df, aes(x = week, y = value), color = "red") +
    xlab("Week of gestation") +
    ylab("Week-specific exposure effect") +
    ylim(-.41, .41)

  lag_plt_all_ns

  # save plot
  ggsave(
    paste0(
      "~/ZeroFillDLM/Simulation/FinalTablesFigures/Lag_figure_mingest",
      min_gest_age,
      ".png"
    ),
    plot = lag_plt_all_ns,
    device = "png",
    width = 7.5,
    height = 7.5,
    dpi = 320
  )

  #-------------------------------------------
  # coverage plots by lag week
  #-------------------------------------------

  # plots by lag week of coverage
  lag_plt_data_all_cover <-
    combined_results_all |>
    select(-starts_with("est")) |>
    filter(model == "ns (best AIC)") |>
    pivot_longer(
      cols = paste0("cover", 1:42),
      names_to = "week",
      values_to = "cover"
    ) |>
    mutate(week = as.numeric(substr(week, 6, 7))) |>
    select(-c("df", "AIC", "min_aic")) |>
    group_by(model, data, GAcontrol, scenario, week) |>
    summarize(
      mean = mean(cover, na.rm = TRUE)
    )
  # this will produce warnings because of the 37 week group with all NAs after week 37

  # remove extra weeks from 37 week truncation so plot doesn't
  # look like the estimates goes to 0 in week 38
  lag_plt_data_all_cover2 <- lag_plt_data_all_cover |>
    filter(
      data != paste0("truncated at week ", min_gest_age) |
        week <= min_gest_age
    )

  # # plots by lag week of coverage with ns only for publication main text
  # lag_plt_all_cover_ns <-
  #   ggplot(
  #     data = lag_plt_data_all_cover2 |>
  #       filter(model %in% c("ns (best AIC)"))
  #   ) +
  #   geom_line(aes(x = week, y = mean)) +
  #   # facet_grid(scenario+GAcontrol ~ data) +
  #   # facet_grid2(
  #   #   rows = vars(scenario, GAcontrol),
  #   #   cols = vars(data),
  #   #   strip = strip_nested()
  #   # ) +
  #   theme_bw(base_size = 9) +
  #   theme(strip.background = element_blank()) +
  #   facet_nested(scenario + GAcontrol ~ data, nest_line = TRUE) +
  #   geom_hline(yintercept = 0.95, linetype = 2) +
  #   xlab("Week of gestation") +
  #   ylab("Coverage")
  #
  # # save plot
  # ggsave(
  #   paste0(
  #     "~/ZeroFillDLM/Simulation/FinalTablesFigures/Lag_coverage_figure_mingest",
  #     min_gest_age,
  #     ".png"
  #   ),
  #   plot = lag_plt_all_cover_ns,
  #   device = "png",
  #   width = 7.5,
  #   height = 7.5,
  #   dpi = 320
  # )

  #-------------------------------------------
  # coverage plots by lag week overlay
  #-------------------------------------------

  lag_plt_all_cover_overlay <-
    ggplot(
      data = lag_plt_data_all_cover2 |>
        filter(model %in% c("ns (best AIC)"))
    ) +
    geom_line(aes(x = week, y = mean, color = data, linetype = data)) +
    # facet_grid(scenario+GAcontrol ~ data) +
    # facet_grid2(
    #   rows = vars(scenario, GAcontrol),
    #   cols = vars(data),
    #   strip = strip_nested()
    # ) +
    theme_bw(base_size = 9) +
    theme(strip.background = element_blank(), legend.position = "bottom") +
    facet_nested(GAcontrol ~ scenario, nest_line = TRUE) +
    geom_hline(yintercept = 0.95, linetype = 2) +
    xlab("Week of gestation") +
    ylab("95% Confidence Interval Coverage") +
    labs(color = "Data method", linetype = "Data method")
  lag_plt_all_cover_overlay

  # save plot
  ggsave(
    paste0(
      "~/ZeroFillDLM/Simulation/FinalTablesFigures/Lag_coverage_figure_overlay_mingest",
      min_gest_age,
      ".png"
    ),
    plot = lag_plt_all_cover_overlay,
    device = "png",
    width = 7.5,
    height = 4,
    dpi = 320
  )
}
