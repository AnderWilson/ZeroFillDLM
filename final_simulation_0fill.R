# Simulation study to assess different methods of data processing
# for a distributed lag analysis with different length gestational
# ages including 0-filling after birth and using post-birth exposures
#
# Written by Ander Wilson
# Written 10/28/2025
# Last Modified 12/9/2025
#
# To run this script on Magpie command line use the following command
# from the folder "ZeroFillDLM/"
# Rscript final_simulation_0fill.R > final_simulation_0fill.out 2>&1 &

rm(list = ls())
gc()
library(splines)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)


# ---------------------------------------- #
# ---------------------------------------- #
# simulation settings
# ---------------------------------------- #
# ---------------------------------------- #

# number of simulated data sets
nsims <- 2000

# sample size for each simulated data set
n <- 1000

# range of df to test for natural splines
df_range <- 3:8

# true exposure effect
truth <- rep(0, 42)
truth[15:30] <- c(1 / (1 + exp(seq(3, -4))), 1 / (1 + exp(seq(-4, 3))))
truth <- -truth / sum(truth) # rescale to have cumulative effect of -1
# plot(truth)

# this is used for plotting later
truth_df <- data.frame(week = 1:42, value = truth)

# load exposure data
exp_dta <- read_csv("Data/ExposureDataPM25.csv")


# ---------------------------------------- #
# ---------------------------------------- #
# functions to summarize results from models
# ---------------------------------------- #
# ---------------------------------------- #

# function to summarize results for natural spline model
summary_ns <- function(fit, B, truth) {
  # point estimate
  est <- B %*% fit$coefficients[c(2:(df + 1))]

  # test if interval will cover 0 for each week
  cover <- (abs(est - truth[1:nrow(B)]) /
    sqrt(diag(
      B %*%
        vcov(fit)[c(2:(df + 1)), c(2:(df + 1))] %*%
        t(B)
    ))) <
    qt(0.975, df = fit$df.residual)

  # test if interval will cover 0 for cumulative effect
  cover_CE <- (abs(sum(est) - sum(truth[1:nrow(B)])) /
    sqrt(sum(
      B %*%
        vcov(fit)[c(2:(df + 1)), c(2:(df + 1))] %*%
        t(B)
    ))) <
    qt(0.975, df = fit$df.residual)

  # compile output
  out <- c(
    AIC(fit), # AIC
    mean(cover), # Coverage for lag function from 1 to end
    ifelse(nrow(B) == 42, mean(cover[38:42]), NA), # coverage for lag function from weeks 38-42
    sum(est), # cumulative effect estimate
    cover_CE, # coverage for cumulative effect
    est # week-specific association
  )

  # add 0s to the 37 week results for point estimate
  if (nrow(B) == 37) {
    out <- c(out, rep(0, 5))
  }

  # add week-specific coverage
  out <- c(out, cover)

  # add 0s to the 37 week results for coverage
  if (nrow(B) == 37) {
    out <- c(out, rep(0, 5))
  }

  return(out)
}


# ---------------------------------------- #
# ---------------------------------------- #
# Run simulation
# ---------------------------------------- #
# ---------------------------------------- #

# loop over scenarios
for (scenario in c("B")) {
  print("------------------------")
  print("------------------------")
  print(paste("Scenario", scenario))
  print("------------------------")
  print("------------------------")
  # A: gest. age indep. of exposure, gest. age indep. of BW
  # B: gest. age indep. of exposure, gest. age effects BW. This is the only scenario used in the manuscript.
  # C: gest. age effects exposure, gest. age effects BW

  # loop over GA control or not
  for (GAcontrol in c(TRUE, FALSE)) {
    # place to store results
    results_full <- results_37 <- results_0 <- results_mean <-
      results_carryforward <-
        results_mean <- cbind(
          data.frame(
            simnum = rep(1:nsims, each = length(df_range)),
            model = rep("ns", length(df_range) * nsims),
            df = rep(df_range, nsims),
            AIC = NA,
            coverage_142 = NA,
            coverage_3842 = NA,
            cumulative = NA,
            cumulative_coverage = NA
          ),
          matrix(0, nsims * length(df_range), 42), # for point estimate
          matrix(0, nsims * length(df_range), 42) # cor coverage
        )

    # loop over replicate data sets
    for (sim in 1:nsims) {
      # set a seed for each data set
      set.seed(158 * sim + 21)

      # ---------------------------------------- #
      # Simulate data
      # ---------------------------------------- #

      # select n exposure histories
      exposure <- as.matrix(sample_n(exp_dta, n))

      # generate gestational ages
      if (scenario %in% c("A", "B")) {
        # gestational age independent of exposure
        gestage <- floor(runif(n, 37, 43))
      } else if (scenario %in% c("C")) {
        # gestational age depends on exposure
        gestage <- -scale(rowMeans(exposure[, 20:30])) + rnorm(n) / 2
        gestage <- apply(gestage, 1, function(x) {
          min(which(quantile(gestage, seq(0, 1, length = 6 + 1))[-1] >= x))
        }) +
          36
      }

      # create imputed histories
      exposure0 <- exposuremean <-
        exposurecarryforward <- exposure
      for (i in 1:n) {
        if (gestage[i] < 42) {
          exposure0[i, (gestage[i] + 1):42] <- 0
          exposuremean[i, (gestage[i] + 1):42] <- mean(exposuremean[
            i,
            -c((gestage[i] + 1):42)
          ])
          exposurecarryforward[i, (gestage[i] + 1):42] <- exposurecarryforward[
            i,
            gestage[i]
          ]
        }
      }

      # create truncated histories
      exposure37 <- exposure[, 1:37]

      # simulate outcome with only a main effect of exposure
      set.seed(158 * sim + 21) # reset seed to get same residuals on all scenarios
      y <-
        exposure %*%
        truth + # exposure effect on outcome (negative)
        5 * rnorm(n)

      # add gestational age main effect for scenarios that include that
      if (scenario %in% c("B", "C")) {
        y <- y + (gestage - 37) / 2 # effect of gestational age on outcome (positive)
      }

      # ---------------------------------------- #
      # fit fixed rank natural spline models
      # using the splines package
      # ---------------------------------------- #

      # fit spline model for multiple df
      # select best df for all data sets at the end of the simulation
      for (df in df_range) {
        # spline basis for DLM
        Bfull <- ns(1:42, intercept = TRUE, df = df)
        B37 <- ns(1:37, intercept = TRUE, df = df)

        # make design matrix
        XBfull <- exposure %*% Bfull
        XB0 <- exposure0 %*% Bfull
        XBmean <- exposuremean %*% Bfull
        XB37 <- exposure37 %*% B37
        XBcarryforward <- exposurecarryforward %*% Bfull

        # estimate model
        if (GAcontrol) {
          fit_full <- lm(y ~ XBfull + poly(gestage, 2))
          fit_0 <- lm(y ~ XB0 + poly(gestage, 2))
          fit_mean <- lm(y ~ XBmean + poly(gestage, 2))
          fit_37 <- lm(y ~ XB37 + poly(gestage, 2))
          fit_carryforward <- lm(y ~ XBcarryforward + poly(gestage, 2))
        } else {
          fit_full <- lm(y ~ XBfull)
          fit_0 <- lm(y ~ XB0)
          fit_mean <- lm(y ~ XBmean)
          fit_37 <- lm(y ~ XB37)
          fit_carryforward <- lm(y ~ XBcarryforward)
        }

        # store results
        row_to_save <- (sim - 1) *
          (length(df_range)) +
          which(df_range == df)

        results_full[row_to_save, -c(1:3)] <- summary_ns(fit_full, Bfull, truth)
        results_0[row_to_save, -c(1:3)] <- summary_ns(fit_0, Bfull, truth)
        results_mean[row_to_save, -c(1:3)] <- summary_ns(fit_mean, Bfull, truth)
        results_37[row_to_save, -c(1:3)] <- summary_ns(fit_37, B37, truth)
        results_carryforward[row_to_save, -c(1:3)] <- summary_ns(
          fit_carryforward,
          Bfull,
          truth
        )
      } # end loop over df
    } # end loop over replicate data sets

    # ---------------------------------------- #
    # combine results and summarize
    # ---------------------------------------- #

    # add column names to differentiate estimates and coverage
    colnames(results_full)[-c(1:8)] <- paste0(
      c(rep("est", 42), rep("cover", 42)),
      c(1:42, 1:42)
    )
    colnames(results_0)[-c(1:8)] <- paste0(
      c(rep("est", 42), rep("cover", 42)),
      c(1:42, 1:42)
    )
    colnames(results_mean)[-c(1:8)] <- paste0(
      c(rep("est", 42), rep("cover", 42)),
      c(1:42, 1:42)
    )
    colnames(results_37)[-c(1:8)] <- paste0(
      c(rep("est", 42), rep("cover", 42)),
      c(1:42, 1:42)
    )
    colnames(results_carryforward)[-c(1:8)] <- paste0(
      c(rep("est", 42), rep("cover", 42)),
      c(1:42, 1:42)
    )

    # combine results from all models
    results_all <- bind_rows(
      results_full |> mutate(data = "natural through week 42"),
      results_0 |> mutate(data = "zero-filled"),
      results_mean |> mutate(data = "mean-filled"),
      results_37 |> mutate(data = "truncated at week 37"),
      results_carryforward |> mutate(data = "carry forward")
    )

    # find best models based on AIC
    results_all <-
      bind_rows(
        results_all,
        results_all |>
          filter(model == "ns") |>
          group_by(simnum, data) |>
          mutate(min_aic = min(AIC)) |>
          filter(AIC == min_aic) |>
          mutate(model = "ns (best AIC)")
      )

    # save results
    save(
      results_all,
      file = paste0(
        "FinalOutput/sim_results_revised_scenario_",
        scenario,
        "_GAcontrol",
        GAcontrol,
        ".rda"
      )
    )
    # load(paste0(
    #   "FinalOutput/sim_results_revised_scenario_",
    #   scenario,
    #   "_GAcontrol",
    #   GAcontrol,
    #   ".rda"
    # ))

    # ---------------------------------------- #
    # summarize cumulative effect (CE) bias results
    # ---------------------------------------- #

    # calculate bias in cumulative effect
    CE_all <- results_all |>
      select(simnum, model, df, data, cumulative, cumulative_coverage) |>
      mutate(cumulative_bias = cumulative - sum(truth))

    # make table
    CE_table <-
      CE_all |>
      filter(model == "ns (best AIC)") |>
      group_by(model, data) |>
      summarise(
        CE_bias = mean(cumulative_bias),
        CE_RMSE = mean(cumulative_bias^2),
        CE_coverage = mean(cumulative_coverage),
        CE_SE = sd(cumulative_bias) / sqrt(nsims),
        CE_t = CE_bias / CE_SE,
        CE_p = 2 * pt(abs(CE_t), df = nsims - 1, lower.tail = FALSE)
      ) |>
      select(model, data, CE_bias, CE_RMSE, CE_coverage, CE_SE, CE_t, CE_p) |>
      data.frame()

    # save results
    write_csv(
      CE_table,
      file = paste0(
        "FinalOutput/cumulative_bias_scenario_",
        scenario,
        "_GAcontrol",
        GAcontrol,
        ".csv"
      )
    )

    # lag specific table
    lag_rmse_table <-
      results_all |>
      select(model, data, simnum, coverage_142, coverage_3842)

    # calculate RMSE
    lag_rmse_table$rmse3842 <-
      sqrt(rowMeans(results_all[, paste0("est", 38:42)]^2))
    lag_rmse_table$rmse142 <-
      sqrt(colMeans((t(results_all[, paste0("est", 1:42)]) - truth)^2))

    # summarize results
    lag_rmse_table <-
      lag_rmse_table |>
      filter(model == "ns (best AIC)") |>
      group_by(model, data) |>
      summarise(
        rmse142 = mean(rmse142),
        rmse3842 = mean(rmse3842),
        cover142 = mean(coverage_142),
        cover3842 = mean(coverage_3842)
      ) |>
      select(model, data, rmse142, rmse3842, cover142, cover3842) |>
      mutate(rmse3842 = ifelse(data == "truncated at week 37", NA, rmse3842)) |>
      data.frame()

    # save results
    write_csv(
      lag_rmse_table,
      file = paste0(
        "FinalOutput/lag_rmse_scenario_",
        scenario,
        "_GAcontrol",
        GAcontrol,
        ".csv"
      )
    )

    # plots by lag week
    lag_plt_data_all <-
      results_all |>
      select(-starts_with("cover")) |>
      filter(model == "ns (best AIC)") |>
      pivot_longer(
        cols = paste0("est", 1:42),
        names_to = "week",
        values_to = "est"
      ) |>
      mutate(week = as.numeric(substr(week, 4, 5))) |>
      select(-c("df", "AIC", "min_aic")) |>
      group_by(model, data, week) |>
      summarize(
        mean = mean(est, na.rm = TRUE),
        min = min(est, na.rm = TRUE),
        max = max(est, na.rm = TRUE),
        q95 = quantile(est, 0.95, na.rm = TRUE),
        q05 = quantile(est, 0.05, na.rm = TRUE)
      )
    # this will produce warnings because of the 37 week group with all NAs after week 37

    lag_rmse <-
      results_all |>
      filter(model == "ns (best AIC)") |>
      select(-c("df", "AIC", "min_aic"))

    # remove extra weeks from 37 week truncation so plot doesn't
    # look like the estimates goes to 0 in week 38
    lag_plt_data_all2 <- lag_plt_data_all |>
      filter(data != "truncated at week 37" | week <= 37)

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
      facet_grid(. ~ data) +
      geom_line(data = truth_df, aes(x = week, y = value), color = "red") +
      theme_bw(base_size = 9) +
      xlab("Week of gestation") +
      ylab("Exposure effect")

    # lag_plt_all_ns

    # save plot
    ggsave(
      paste0(
        "FinalOutput/lag_plots_scenario_nsonly_",
        scenario,
        "_GAcontrol",
        GAcontrol,
        ".pdf"
      ),
      plot = lag_plt_all_ns,
      width = 6.5,
      height = 2,
      units = "in"
    )

    # plots by lag week of coverage
    lag_plt_data_all_cover <-
      results_all |>
      select(-starts_with("est")) |>
      filter(model == "ns (best AIC)") |>
      pivot_longer(
        cols = paste0("cover", 1:42),
        names_to = "week",
        values_to = "cover"
      ) |>
      mutate(week = as.numeric(substr(week, 6, 7))) |>
      select(-c("df", "AIC", "min_aic")) |>
      group_by(model, data, week) |>
      summarize(
        mean = mean(cover, na.rm = TRUE)
      )
    # this will produce warnings because of the 37 week group with all NAs after week 37

    # remove extra weeks from 37 week truncation so plot doesn't
    # look like the estimates goes to 0 in week 38
    lag_plt_data_all_cover2 <- lag_plt_data_all_cover |>
      filter(data != "truncated at week 37" | week <= 37)

    # plots by lag week of coverage with ns only for publication main text
    lag_plt_all_cover_ns <-
      ggplot(
        data = lag_plt_data_all_cover2 |> filter(model %in% c("ns (best AIC)"))
      ) +
      geom_line(aes(x = week, y = mean)) +
      facet_grid(. ~ data) +
      ylim(0, 1) +
      theme_bw(base_size = 9) +
      geom_hline(yintercept = 0.95, linetype = 2) +
      xlab("Week of gestation") +
      ylab("Coverage")

    # lag_plt_all

    # save plot
    ggsave(
      paste0(
        "FinalOutput/lag_plots_coverage_scenario_nsonly_",
        scenario,
        "_GAcontrol",
        GAcontrol,
        ".pdf"
      ),
      plot = lag_plt_all_cover_ns,
      width = 6.5,
      height = 2,
      units = "in"
    )

    # make table summarizing both cumulative and lag-specific estiamtes
    final_table <-
      left_join(lag_rmse_table, CE_table, by = c("model", "data")) |>
      mutate(
        CE_bias = CE_bias * 100,
        CE_RMSE = CE_RMSE * 100,
        CE_coverage = CE_coverage * 100,
        RMSE_1_42 = rmse142 * 100,
        RMSE_38_42 = rmse3842 * 100,
        cover_1_42 = cover142 * 100,
        cover_38_42 = cover3842 * 100
      ) |>
      select(
        model,
        data,
        CE_bias,
        CE_RMSE,
        CE_coverage,
        CE_p,
        RMSE_1_42,
        RMSE_38_42,
        cover_1_42,
        cover_38_42
      )

    # save table
    write_csv(
      final_table,
      paste0(
        "FinalOutput/final_table_scenario",
        scenario,
        "_GAcontrol",
        GAcontrol,
        ".csv"
      ),
    )
  } # end loop over GAcontrol
} # end loop over scenario


# ---------------------------------------- #
# ---------------------------------------- #
# create combined final table that includes
# all scenarios
# ---------------------------------------- #
# ---------------------------------------- #

# find the files containing simulation results
filenames <- list.files(path = "FinalOutput/")
filenames <- filenames[grep("final_table_scenario", filenames)]

# load and combine the tables
combined_final_table <- NULL
for (i in filenames) {
  tabletemp <- read_csv(paste0("FinalOutput/", i))
  tabletemp$file <- i
  combined_final_table <- bind_rows(combined_final_table, tabletemp)
}

# save results
write_csv(
  combined_final_table,
  "FinalOutput/combined_final_table.csv",
)
