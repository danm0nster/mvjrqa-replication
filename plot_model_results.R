## ---------------------------
##
## Script name: plot_model_results.R
##
## Purpose of script: Plot results from simulation models.
##
## Author: Dan Moenster & Sebastian Wallot
##
## Date Created: May 26 2025
##
## ---------------------------
##
## Notes:
##
## This script will produce plots used in figures in the paper.
## It uses data produced by the data generation scripts generate_*_data.R.
##
##
## ---------------------------

library(readr)
library(ggplot2)
source("R/utils.R")
source("R/bootstrap.R")
source("R/plot_utils.R")

# Make sure "./plots" directory exists
create_dir_if_not_present("plots")

#
# Linear stochastic system
#

linear_model_coupling_sweep <- read_csv(
  "data/linear_stochastic_coupling_sweep.csv"
)

lin_jrr_rr_plot <- sync_plot(linear_model_coupling_sweep) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  xlab("coupling strength k")
plot(lin_jrr_rr_plot)

ggsave(plot = lin_jrr_rr_plot,
       filename = "plots/linear_stochastic_jrr_rr.pdf",
       width = 1.5 * 4, height = 1.5 * 3)


#
# Logistic maps with external driver
#

logistic_model_coupling_sweep <- read_csv(
  "data/driven_coupled_logistic_maps_coupling_sweep.csv"
)

# TODO: Find values that are not too close, so curves do not clutter
# Maybe 4, 8, 16?
log_jrr_rr_plot_data <-  logistic_model_coupling_sweep |>
  filter(RR_target %in% c(4, 8, 16))


log_jrr_rr_plot <- sync_plot(log_jrr_rr_plot_data) +
  scale_x_continuous(breaks = c(0, 0.02, 0.04, 0.06, 0.08)) +
  xlab(expression(coupling ~ strength ~ eta))

plot(log_jrr_rr_plot)

ggsave(plot = log_jrr_rr_plot,
       filename = "plots/logistic_driven_jrr_rr.pdf",
       width = 1.5 * 4, height = 1.5 * 3)


#
# Lorenz - harmonic oscillator system
#

lorenz_harmonic_coupling_sweep <- read_csv(
  "data/lorenz_harmonic_coupling_sweep.csv"
)

lor_jrr_rr_plot <- sync_plot(lorenz_harmonic_coupling_sweep)
plot(lor_jrr_rr_plot)

ggsave(plot = lor_jrr_rr_plot,
       filename = "plots/lorenz_harmonic_jrr_by_rr.pdf",
       width = 1.5 * 4, height = 1.5 * 3)


#
# Lorenz 96 system coupled to harmonic oscillator
#

lorenz96_coupling_sweep <- read_csv(
  "data/lorenz96_harmonic_coupling_sweep.csv"
)

lor96_jrr_rr_plot <- sync_plot(lorenz96_coupling_sweep) +
  xlab(expression(coupling ~ strength ~ kappa))
plot(lor96_jrr_rr_plot)

ggsave(plot = lor96_jrr_rr_plot,
       filename = "plots/lorenz96_harmonic_jrr_rr.pdf",
       width = 1.5 * 4, height = 1.5 * 3)


#
# The Joint Recurrence Coupling plots
#

#
# Linear stochastic system
#

linear_model_rr_sweep <- read_csv("data/linear_stochastic_rr_sweep.csv")

linear_coupling_colours <- c(
  "0" = "#fd8d3c",
  "3" = "#fc4e2a",
  "6" = "#e31a1c",
  "12" = "#b10026"
)
lin_jrc_plot <- jrc_plot(linear_model_rr_sweep,
                         errorbars = "CI",
                         ci_level = 0.95,
                         coupling_name = "k") +
  scale_color_manual(values = linear_coupling_colours,
                     aesthetics = c("colour", "fill")) +
  scale_x_log10(breaks = c(0.1, 1, 10, 100),
                labels = c(0.1, 1, 10, 100)) +
  scale_y_continuous(breaks = c(0.00, 0.03, 0.06, 0.09, 0.12),
                     labels = c(0.00, 0.03, 0.06, 0.09, 0.12)) +
  theme(text = element_text(size = 22))
plot(lin_jrc_plot)
ggsave(lin_jrc_plot,
       filename = "plots/linear_jrr_by_rr2_log_color.pdf",
       width = 1.5 * 4, height = 1.5 * 3)

#
# Driven logistic map
#

map_coupling_colours <- c(
  "0" = "#fd8d3c",
  "0.02" = "#fc4e2a",
  "0.04" = "#e31a1c",
  "0.08" = "#b10026"
)

logistic_model_data <- read_csv(
  "data/driven_coupled_logistic_maps_rr_sweep.csv"
)

# Calculate the minimum recurrence rate for the periodic driver
# and add a vertical line at the minimum possible RR
# The recurrence plot is 500 x 500
n_rp <- 500
# The driving signal has a period of 30
period <- 30
rec_per_row <- floor((n_rp - 1) / period)
rr_periodic <- n_rp * rec_per_row / n_rp^2
rr_periodic_pct <- 100 * rr_periodic

map_jrc_plot <- jrc_plot(logistic_model_data,
                         errorbars = "CI",
                         ci_level = 0.95,
                         coupling_name = "eta",
                         show_zero = FALSE) +
  scale_color_manual(values = map_coupling_colours,
                     aesthetics = c("colour", "fill")) +
  scale_x_log10(breaks = c(3, 10, 30, 100), labels = c(3, 10, 30, 100)) +
  scale_y_continuous(breaks = c(0.00, 0.03, 0.06, 0.09, 0.12),
                     labels = c(0.00, 0.03, 0.06, 0.09, 0.12)) +
  theme(text = element_text(size = 22)) +
  geom_vline(xintercept = rr_periodic_pct,
             linetype = "dashed",
             color = "grey") +
  expand_limits(y = c(0, 0.12)) # Manual adjustment of scale

plot(map_jrc_plot)

ggsave(plot = map_jrc_plot,
       filename = "plots/logistic_maps_external_driver_color.pdf",
       width = 1.5 * 4, height = 1.5 * 3)

#
# Lorenz & Harmonic oscillator
#

lorenz_coupling_colours <- c(
  "0" = "#fd8d3c",
  "0.1" = "#fc4e2a",
  "0.2" = "#e31a1c",
  "0.4" = "#b10026"
)

lorenz_harmonic_rr_sweep <- read_csv("data/lorenz_harmonic_rr_sweep.csv")

lor_jrc_plot <- jrc_plot(lorenz_harmonic_rr_sweep,
                         errorbars = "CI",
                         ci_level = 0.95,
                         coupling_name = "c") +
  scale_color_manual(values = lorenz_coupling_colours,
                     aesthetics = c("colour", "fill")) +
  expand_limits(x = 0, y = 1.4) +
  scale_x_log10(breaks = c(0.1, 1, 10, 100), labels = c(0.1, 1, 10, 100)) +
  scale_y_continuous(
    breaks = seq(0, 1.25, by = 0.25),
    labels = c("0.00", "0.25", "0.50", "0.75", "1.00", "1.25")
  ) +
  theme(text = element_text(size = 22))


plot(lor_jrc_plot)

ggsave(lor_jrc_plot,
       filename = "plots/lor-osc_jrr_by_rr2_log_color.pdf",
       width = 1.5 * 4, height = 1.5 * 3)


#
# Lorenz 96 system coupled to harmonic oscillator
#

lorenz96_coupling_colours <- c(
  "0" = "#fd8d3c",
  "0.1" = "#fc4e2a",
  "0.2" = "#e31a1c",
  "0.4" = "#b10026"
)

lorenz96_harmonic_rr_sweep <- read_csv(
  "data/lorenz96_harmonic_rr_sweep.csv"
)

lor_96_jrc_data <- lorenz96_harmonic_rr_sweep |>
  filter(coupling %in% c(0, 0.1, 0.2, 0.4))

lor96_jrc_plot <- jrc_plot(lor_96_jrc_data,
                           errorbars = "CI",
                           ci_level = 0.95,
                           coupling_name = "kappa") +
  scale_color_manual(values = lorenz96_coupling_colours,
                     aesthetics = c("colour", "fill")) +
  expand_limits(x = 0, y = 1.4) +
  scale_x_log10(breaks = c(0.1, 1, 10, 100), labels = c(0.1, 1, 10, 100)) +
  scale_y_continuous(
    breaks = seq(0, 1.25, by = 0.25),
    labels = c("0.00", "0.25", "0.50", "0.75", "1.00", "1.25")
  ) +
  theme(text = element_text(size = 22))

plot(lor96_jrc_plot)

ggsave(plot = lor96_jrc_plot,
       filename = "plots/lor96_osc_jrr_by_rr2_log_color.pdf",
       width = 1.5 * 4, height = 1.5 * 3)
