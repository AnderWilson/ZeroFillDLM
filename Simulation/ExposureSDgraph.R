# Make graphs of the exposure standard deviation and correlation
# with gestational age
#
# Written by Ander Wilson
# Written 10/29/2025
# Last Modified 4/23/2026
#
# To run this script on Magpie command line use the following command
# from the folder "ZeroFillDLM/"
# Rscript ExposureSDgraph.R > ExposureSDgraph.out 2>&1 &

#--------------------------------------------------
# load exposure data
#--------------------------------------------------

# use all exposure data
exposure <- read_csv("Data/ExposureDataPM25.csv")
exposure <- as.matrix(exposure)[, 1:42]
n <- nrow(exposure)

# use subsample, which is faster but more variable
# n <- 1000
#
# exposure <- as.matrix(sample_n(exp_dta, n))

#--------------------------------------------------
# generate gestational ages
#--------------------------------------------------

gestdist <- c(
  0.001,
  0.004,
  0.007,
  0.011,
  0.020,
  0.033,
  0.063,
  0.142,
  0.298,
  0.656,
  0.901,
  0.995,
  1.000
)
# above is cumulative, get the week-specific probability
gestdist <- gestdist - c(0, gestdist[-length(gestdist)])

names(gestdist) <- 30:42

gestage <- sample(
  37:42,
  size = n,
  prob = gestdist[as.character(37:42)],
  replace = TRUE
)

table(gestage)

# use this for uniform gestational age, which is unrealistic
# gestage <- floor(runif(n, 37, 43))

#--------------------------------------------------
# create imputed histories
#--------------------------------------------------
exposure0 <- exposuremean <-
  exposurenoise <- exposurecarryforward <- exposure
for (i in 1:n) {
  if (gestage[i] < 42) {
    exposure0[i, (gestage[i] + 1):42] <- 0
    exposuremean[i, (gestage[i] + 1):42] <- mean(exposuremean[
      i,
      -c((gestage[i] + 1):42)
    ])
    exposurenoise[i, (gestage[i] + 1):42] <- sample(
      exposure,
      length((gestage[i] + 1):42),
      replace = TRUE
    )
    exposurecarryforward[i, (gestage[i] + 1):42] <- exposurecarryforward[
      i,
      gestage[i]
    ]
  }
}

# create truncated histories
exposure37 <- exposure[, 1:37]

# plot data for standard deviation
sd_data <-
  data.frame(
    week = c(rep(1:42, 4), 1:37),
    sd = c(
      apply(exposure0, 2, sd),
      apply(exposure, 2, sd),
      apply(exposuremean, 2, sd),
      apply(exposurecarryforward, 2, sd),
      apply(exposure37, 2, sd)
    ),
    data = c(
      rep(
        c("zero-fill", "natural through week 42", "mean-fill", "carry forward"),
        each = 42
      ),
      rep("truncated at week 37", 37)
    )
  )

#--------------------------------------------------
# plot for standard deviation
#--------------------------------------------------
sd_plt <-
  sd_data |>
  ggplot(aes(x = week, y = sd, color = data, linetype = data)) +
  # geom_line(linewidth = .5) + # use for pdf
  # theme_bw(base_size = 7) + # use for pdf
  geom_line(linewidth = 1) + # use for eps
  theme_bw(base_size = 10) + # use for eps
  ylab("Week-specific exposure\nstandard deviation") +
  xlab("Week of gestation") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(.05, .05),
    legend.justification = c(0, 0),
    legend.title = element_blank(),
    legend.box.background = element_rect(colour = "black"),
    # legend.key.size = unit(.12, 'in'), # use for pdf
    # legend.key.width = unit(.3, 'in') # use for pdf
    legend.key.size = unit(.16, 'in'), # use for eps
    legend.key.width = unit(.5, 'in') # use for eps
  ) +
  scale_linetype_manual(values = c(2, 3, 4, 5, 1))
# ggtitle("Exposure Standard Deviation by Week")
# sd_plt

# ggsave(
#   "Simulation/FinalOutput/mingest37/exposure_sd_graph.pdf",
#   plot = sd_plt,
#   width = 3,
#   height = 2,
#   units = "in"
# )

ggsave(
  "~/ZeroFillDLM/Simulation/FinalTablesFigures/exposure_sd_graph.pmg",
  plot = sd_plt,
  device = "png",
  width = 4,
  height = 3,
  units = "in",
  dpi = 320
)


# plot data for correlation with gestational age
cor_data <-
  data.frame(
    week = c(rep(1:42, 4), 1:37),
    cor = c(
      apply(exposure0, 2, function(x) {
        cor(x, gestage)
      }),
      apply(exposure, 2, function(x) {
        cor(x, gestage)
      }),
      apply(exposuremean, 2, function(x) {
        cor(x, gestage)
      }),
      apply(exposurecarryforward, 2, function(x) {
        cor(x, gestage)
      }),
      apply(exposure37, 2, function(x) {
        cor(x, gestage)
      })
    ),
    data = c(
      rep(
        c("zero-fill", "natural through week 42", "mean-fill", "carry forward"),
        each = 42
      ),
      rep("truncated at week 37", 37)
    )
  )

#--------------------------------------------------
# plot for correlation with gestational age
#--------------------------------------------------
cor_plt <-
  cor_data |>
  ggplot(aes(x = week, y = cor, color = data, linetype = data)) +
  # geom_line(linewidth = .5) + # use for pdf
  # theme_bw(base_size = 7) + # use for pdf
  geom_line(linewidth = 1) + # use for eps
  theme_bw(base_size = 10) + # use for eps
  # geom_hline(yintercept = c(0, 1), linewidth=.2) +
  ylab(
    "Correlation between week-specific\nexposure and gestational age at birth"
  ) +
  xlab("Week of gestation") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(.05, .95),
    legend.justification = c(0, 1),
    legend.title = element_blank(),
    legend.box.background = element_rect(colour = "black"),
    # legend.key.size = unit(.12, 'in'), # use for pdf
    # legend.key.width = unit(.3, 'in') # use for pdf
    legend.key.size = unit(.16, 'in'), # use for eps
    legend.key.width = unit(.5, 'in') # use for eps
  ) +
  scale_linetype_manual(values = c(2, 3, 4, 5, 1))
# ggtitle("Correlation between exposure and gestational age by Week")
cor_plt

ggsave(
  "~/ZeroFillDLM/Simulation/FinalTablesFigures/exposure_gestage_cor_graph.png",
  plot = cor_plt,
  device = "png",
  width = 4,
  height = 3,
  units = "in",
  dpi = 320
)
