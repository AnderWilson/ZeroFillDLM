# Make graphs of the exposure standard deviation and correlation
# with gestational age
#
# Written by Ander Wilson
# Written 10/29/2025
# Last Modified 12/9/2025
#
# To run this script on Magpie command line use the following command
# from the folder "ZeroFillDLM/"
# Rscript ExposureSDgraph.R > ExposureSDgraph.out 2>&1 &

# use all exposure data
exposure <- read_csv("Data/ExposureDataPM25.csv")
exposure <- as.matrix(exposure)
n <- nrow(exposure)

# n <- 10000
#
# exposure <- as.matrix(sample_n(exp_dta, n))

# generate gestational ages
gestage <- floor(runif(n, 37, 43))

# create imputed histories
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

# plot for standard deviation
sd_plt <-
  sd_data |>
  ggplot(aes(x = week, y = sd, color = data, linetype = data)) +
  geom_line(linewidth = .5) +
  theme_bw(base_size = 7) +
  # geom_hline(yintercept = 0) +
  # ylab("Stadnard Deviation") +
  ylab("Week-specific exposure\nstandard deviation") +
  xlab("Week of gestation") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(.05, .95),
    legend.justification = c(0, 1),
    legend.title = element_blank(),
    legend.box.background = element_rect(colour = "black"),
    legend.key.size = unit(.12, 'in'),
    legend.key.width = unit(.3, 'in')
  ) +
  scale_linetype_manual( values=c(2,3,4,5,1))
# ggtitle("Exposure Standard Deviation by Week")
sd_plt

ggsave(
  "FinalOutput/exposure_sd_graph.pdf",
  plot = sd_plt,
  width = 3,
  height = 2,
  units = "in"
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

# plot for correlation with gestational age
cor_plt <-
  cor_data |>
  ggplot(aes(x = week, y = cor, color = data, linetype = data)) +
  geom_line(linewidth = .5) +
  theme_bw(base_size = 7) +
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
    legend.key.size = unit(.12, 'in'),
    legend.key.width = unit(.3, 'in')
  ) +
  scale_linetype_manual( values=c(2,3,4,5,1))
# ggtitle("Correlation between exposure and gestational age by Week")
cor_plt

ggsave(
  "FinalOutput/exposure_gestage_cor_graph.pdf",
  plot = cor_plt,
  width = 3,
  height = 2,
  units = "in"
)
