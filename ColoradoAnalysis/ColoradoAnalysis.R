# this script reads in the Colorado birth data
# and strips off only the exposure data
# and gestational age at birth

# from the folder "ZeroFillDLM/ColoradoAnalysis"
# Rscript ColoradoAnalysis.R > ColoradoAnalysis.out 2>&1 &

library(dplyr)
library(tidyr)
library(readr)
library(splines)
library(ggplot2)
library(patchwork)

rm(list = ls())
gc()


# make figure showing distribution of gestational age at birth

load("~/ZeroFillDLM/Data/co_dlm_data.rda")

p_gestage <-
  ggplot(co_birth_dlm_data, aes(x = EstGest)) +
  geom_histogram(
    aes(y = after_stat(count / sum(count))),
    binwidth = 1,
    color = "gray70",
    fill = "gray70"
  ) +
  theme_bw() +
  labs(
    x = "Gestational age at birth",
    y = "Proportion of births (weeks)"
  )


ggsave(
  file = "~/ZeroFillDLM/ColoradoAnalysis/FiguresTables/colorado_gestage.png",
  plot = p_gestage,
  device = "png",
  height = 3,
  width = 5,
  units = "in",
  dpi = 320
)


# range of df to test for natural splines
df_range <- 4

for (mingest in c(30, 37)) {
  # load data
  load("~/ZeroFillDLM/Data/co_dlm_data.rda")

  # filter to only full term births.
  dta <-
    co_birth_dlm_data |>
    filter(EstGest >= mingest)

  # some
  n <- nrow(dta)
  gestage <- dta |> pull(EstGest)

  # get exposure data
  exposure <-
    dta |>
    select(starts_with("cmaq_pm25_")) |>
    as.matrix()

  # simplify column names
  colnames(exposure) <- paste0("pm25_", 1:42)

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
  min(gestage)
  exposure37 <- exposure[, 1:min(gestage)]

  # place to store results
  model_fits <- list(
    full_GA = list(),
    mean_GA = list(),
    zero_GA = list(),
    truncated_GA = list(),
    carryforward_GA = list(),
    full_noGA = list(),
    mean_noGA = list(),
    zero_noGA = list(),
    truncated_noGA = list(),
    carryforward_noGA = list()
  )

  # AIC results (not used in final run)
  AICresults <- matrix(NA, 10, max(df_range))
  rownames(AICresults) <- names(model_fits)

  for (df in df_range) {
    # spline basis for DLM
    Bfull <- ns(1:42, intercept = TRUE, df = df)
    B37 <- ns(1:min(gestage), intercept = TRUE, df = df)

    # make design matrix
    XBfull <- exposure %*% Bfull
    XB0 <- exposure0 %*% Bfull
    XBmean <- exposuremean %*% Bfull
    XB37 <- exposure37 %*% B37
    XBcarryforward <- exposurecarryforward %*% Bfull

    # estimate model
    # models with GA control
    fit_full <- lm(
      BWGr ~
        XBfull +
          poly(EstGest, 2) +
          Sex +
          as.factor(MOC) *
            as.factor(YOC) +
          elev_feet_tract +
          poly(MatAge, 2) +
          Marital2 +
          Income +
          MEduc +
          PrenatalCare +
          MotherBMI +
          hispanic +
          race +
          Smk +
          fipscoor +
          tmmx_tri1 +
          tmmx_tri2 +
          tmmx_tri3,
      data = dta
    )
    AICresults["full_GA", df] <- AIC(fit_full)
    model_fits$full_GA[[df]] <- list(
      fit = fit_full,
      df = df,
      basis = Bfull,
      GA = TRUE,
      name = "natural through week 42"
    )

    fit_0 <- lm(
      BWGr ~
        XB0 +
          poly(EstGest, 2) +
          Sex +
          as.factor(MOC) *
            as.factor(YOC) +
          elev_feet_tract +
          poly(MatAge, 2) +
          Marital2 +
          Income +
          MEduc +
          PrenatalCare +
          MotherBMI +
          hispanic +
          race +
          Smk +
          fipscoor +
          tmmx_tri1 +
          tmmx_tri2 +
          tmmx_tri3,
      data = dta
    )
    AICresults["zero_GA", df] <- AIC(fit_0)
    model_fits$zero_GA[[df]] <- list(
      fit = fit_0,
      df = df,
      basis = Bfull,
      GA = TRUE,
      name = "zero-filled"
    )

    fit_mean <- lm(
      BWGr ~
        XBmean +
          poly(EstGest, 2) +
          Sex +
          as.factor(MOC) *
            as.factor(YOC) +
          elev_feet_tract +
          poly(MatAge, 2) +
          Marital2 +
          Income +
          MEduc +
          PrenatalCare +
          MotherBMI +
          hispanic +
          race +
          Smk +
          fipscoor +
          tmmx_tri1 +
          tmmx_tri2 +
          tmmx_tri3,
      data = dta
    )
    AICresults["mean_GA", df] <- AIC(fit_mean)
    model_fits$mean_GA[[df]] <- list(
      fit = fit_mean,
      df = df,
      basis = Bfull,
      GA = TRUE,
      name = "mean-filled"
    )

    fit_37 <- lm(
      BWGr ~
        XB37 +
          poly(EstGest, 2) +
          Sex +
          as.factor(MOC) *
            as.factor(YOC) +
          elev_feet_tract +
          poly(MatAge, 2) +
          Marital2 +
          Income +
          MEduc +
          PrenatalCare +
          MotherBMI +
          hispanic +
          race +
          Smk +
          fipscoor +
          tmmx_tri1 +
          tmmx_tri2 +
          tmmx_tri3,
      data = dta
    )
    AICresults["truncated_GA", df] <- AIC(fit_37)
    model_fits$truncated_GA[[df]] <- list(
      fit = fit_37,
      df = df,
      basis = B37,
      GA = TRUE,
      name = "truncated at week 37"
    )

    fit_carryforward <- lm(
      BWGr ~
        XBcarryforward +
          poly(EstGest, 2) +
          Sex +
          as.factor(MOC) *
            as.factor(YOC) +
          elev_feet_tract +
          poly(MatAge, 2) +
          Marital2 +
          Income +
          MEduc +
          PrenatalCare +
          MotherBMI +
          hispanic +
          race +
          Smk +
          fipscoor +
          tmmx_tri1 +
          tmmx_tri2 +
          tmmx_tri3,
      data = dta
    )
    AICresults["carryforward_GA", df] <- AIC(fit_carryforward)
    model_fits$carryforward_GA[[df]] <- list(
      fit = fit_carryforward,
      df = df,
      basis = Bfull,
      GA = TRUE,
      name = "carry forward"
    )

    #models with no GA control

    fit_full <- lm(
      BWGr ~
        XBfull +
          Sex +
          as.factor(MOC) *
            as.factor(YOC) +
          elev_feet_tract +
          poly(MatAge, 2) +
          Marital2 +
          Income +
          MEduc +
          PrenatalCare +
          MotherBMI +
          hispanic +
          race +
          Smk +
          fipscoor +
          tmmx_tri1 +
          tmmx_tri2 +
          tmmx_tri3,
      data = dta
    )
    AICresults["full_noGA", df] <- AIC(fit_full)
    model_fits$full_noGA[[df]] <- list(
      fit = fit_full,
      df = df,
      basis = Bfull,
      GA = FALSE,
      name = "natural through week 42"
    )

    fit_0 <- lm(
      BWGr ~
        XB0 +
          Sex +
          as.factor(MOC) *
            as.factor(YOC) +
          elev_feet_tract +
          poly(MatAge, 2) +
          Marital2 +
          Income +
          MEduc +
          PrenatalCare +
          MotherBMI +
          hispanic +
          race +
          Smk +
          fipscoor +
          tmmx_tri1 +
          tmmx_tri2 +
          tmmx_tri3,
      data = dta
    )
    AICresults["zero_noGA", df] <- AIC(fit_0)
    model_fits$zero_noGA[[df]] <- list(
      fit = fit_0,
      df = df,
      basis = Bfull,
      GA = FALSE,
      name = "zero-filled"
    )

    fit_mean <- lm(
      BWGr ~
        XBmean +
          Sex +
          as.factor(MOC) *
            as.factor(YOC) +
          elev_feet_tract +
          poly(MatAge, 2) +
          Marital2 +
          Income +
          MEduc +
          PrenatalCare +
          MotherBMI +
          hispanic +
          race +
          Smk +
          fipscoor +
          tmmx_tri1 +
          tmmx_tri2 +
          tmmx_tri3,
      data = dta
    )
    AICresults["mean_noGA", df] <- AIC(fit_mean)
    model_fits$mean_noGA[[df]] <- list(
      fit = fit_mean,
      df = df,
      basis = Bfull,
      GA = FALSE,
      name = "mean-filled"
    )

    fit_37 <- lm(
      BWGr ~
        XB37 +
          Sex +
          as.factor(MOC) *
            as.factor(YOC) +
          elev_feet_tract +
          poly(MatAge, 2) +
          Marital2 +
          Income +
          MEduc +
          PrenatalCare +
          MotherBMI +
          hispanic +
          race +
          Smk +
          fipscoor +
          tmmx_tri1 +
          tmmx_tri2 +
          tmmx_tri3,
      data = dta
    )
    AICresults["truncated_noGA", df] <- AIC(fit_37)
    model_fits$truncated_noGA[[df]] <- list(
      fit = fit_37,
      df = df,
      basis = B37,
      GA = FALSE,
      name = "truncated at week 37"
    )

    fit_carryforward <- lm(
      BWGr ~
        XBcarryforward +
          Sex +
          as.factor(MOC) *
            as.factor(YOC) +
          elev_feet_tract +
          poly(MatAge, 2) +
          Marital2 +
          Income +
          MEduc +
          PrenatalCare +
          MotherBMI +
          hispanic +
          race +
          Smk +
          fipscoor +
          tmmx_tri1 +
          tmmx_tri2 +
          tmmx_tri3,
      data = dta
    )
    AICresults["carryforward_noGA", df] <- AIC(fit_carryforward)
    model_fits$carryforward_noGA[[df]] <- list(
      fit = fit_carryforward,
      df = df,
      basis = Bfull,
      GA = FALSE,
      name = "carry forward"
    )
  } # end loop over df

  #--------------------------------------------------
  # make table with cumulative effect estimates
  #--------------------------------------------------

  plot_data <- NULL
  cumulative_effects <- NULL

  for (i in rownames(AICresults)) {
    df_min <- which.min(AICresults[i, ])

    df <- model_fits[[i]][[df_min]]$df
    basis <- model_fits[[i]][[df_min]]$basis
    beta <- coef(model_fits[[i]][[df_min]]$fit)[2:(df + 1)]
    beta_var <- vcov(model_fits[[i]][[df_min]]$fit)[2:(df + 1), 2:(df + 1)]
    results_temp <-
      data.frame(
        model = i,
        model_name = model_fits[[i]][[df_min]]$name,
        gest_control_logical = model_fits[[i]][[df_min]]$GA,
        gest_control = ifelse(
          model_fits[[i]][[df_min]]$GA,
          "adjusted for gestational age",
          "not adjusted for gestational age"
        ),
        week = 1:nrow(basis),
        estimate = basis %*% beta,
        se = sqrt(diag(basis %*% beta_var %*% t(basis)))
      )
    results_temp$lower <- results_temp$estimate -
      results_temp$se * qt(0.975, model_fits[[i]][[df_min]]$fit$df.residual)
    results_temp$upper <- results_temp$estimate +
      results_temp$se * qt(0.975, model_fits[[i]][[df_min]]$fit$df.residual)

    plot_data <- rbind(plot_data, results_temp)

    cumulative_effects_temp <-
      data.frame(
        model = i,
        model_name = model_fits[[i]][[df_min]]$name,
        gest_control_logical = model_fits[[i]][[df_min]]$GA,
        gest_control = ifelse(
          model_fits[[i]][[df_min]]$GA,
          "adjusted for gestational age",
          "not adjusted for gestational age"
        ),
        week = "cumulative",
        estimate = sum(basis %*% beta),
        se = sqrt(sum(basis %*% beta_var %*% t(basis)))
      )
    cumulative_effects_temp$lower <- cumulative_effects_temp$estimate -
      cumulative_effects_temp$se *
        qt(0.975, model_fits[[i]][[df_min]]$fit$df.residual)
    cumulative_effects_temp$upper <- cumulative_effects_temp$estimate +
      cumulative_effects_temp$se *
        qt(0.975, model_fits[[i]][[df_min]]$fit$df.residual)

    cumulative_effects <- rbind(cumulative_effects, cumulative_effects_temp)
  }

  cumulative_1 <- cumulative_effects |>
    filter(gest_control_logical) |>
    select(
      model_name,
      GA_estimate = estimate,
      GA_se = se,
      GA_lower = lower,
      GA_upper = upper
    )
  cumulative_2 <- cumulative_effects |>
    filter(!gest_control_logical) |>
    select(
      model_name,
      noGA_estimate = estimate,
      noGA_se = se,
      noGA_lower = lower,
      noGA_upper = upper
    )
  cumulative_table <- full_join(cumulative_1, cumulative_2, by = "model_name")
  
  cumulative_table <- cumulative_table |> arrange(model_name)

  # table with cumulative effect estimates
  write_csv(
    cumulative_table,
    paste0(
      "~/ZeroFillDLM/ColoradoAnalysis/FiguresTables/colorado_analysis_cumulative_association_table",
      mingest,
      ".csv"
    )
  )

  #--------------------------------------------------
  # make figure with all exposure-time-response functions
  #--------------------------------------------------

  plot_panels <- list()
  letter <- 1
  for (j in unique(plot_data$gest_control)[order(unique(
    plot_data$gest_control
  ))]) {
    for (i in unique(plot_data$model_name)[order(unique(
      plot_data$model_name
    ))]) {
      tempdata <- plot_data |>
        filter(model_name == i & gest_control == j)
      lwr <- min(-1, min(tempdata$lower) - 0.1)
      upr <- max(.5, max(tempdata$upper) + 0.1)

      plot_panels[[paste(i, j)]] <-
        ggplot(
          data = tempdata,
          aes(x = week, y = estimate, ymin = lower, ymax = upper)
        ) +
        geom_hline(yintercept = 0, linetype = 3) +
        geom_ribbon(fill = "grey80", alpha = .7) +
        geom_line() +
        theme_bw(base_size = 9) +
        theme(plot.title = element_text(size = 6)) +
        ylim(lwr, upr) +
        xlim(0, 42) +
        ylab(expression(atop(
          "Estimated difference in birth weight (g)",
          "per 1 " * mu * "g/" * m^3 * " increase in PM"[2.5] * " exposures"
        ))) +
        xlab("Week of gestation") +
        ggtitle(paste0(LETTERS[letter], ") ", i, ",\n", j))

      letter <- letter + 1
    }
  }

  p <- plot_panels[[1]] +
    plot_panels[[2]] +
    plot_panels[[3]] +
    plot_panels[[4]] +
    plot_panels[[5]] +
    plot_panels[[6]] +
    plot_panels[[7]] +
    plot_panels[[8]] +
    plot_panels[[9]] +
    plot_panels[[10]] +
    plot_layout(axes = "collect", ncol = 5)

  ggsave(
    paste0(
      "~/ZeroFillDLM/ColoradoAnalysis/FiguresTables/colorado_full_analysis_mingest",
      mingest,
      ".png"
    ),
    plot = p,
    device = "png",
    height = 4,
    width = 8,
    units = "in",
    dpi = 320
  )

  # sample size
  nrow(dta)
  write(
    paste("n=", nrow(dta)),
    paste0(
      "~/ZeroFillDLM/ColoradoAnalysis/FiguresTables/colorado_analysis_sample_size_mingest",
      mingest,
      ".txt"
    )
  )
}
