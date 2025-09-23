## ---------------------------
##
## Script name: plot_extreme_coupling.R
##
## Purpose of script: Plot JRR/RR for extreme values of coupling.
##
## Author: Dan Moenster & Sebastian Wallot
##
## Date Created: May 27 2025
##
## ---------------------------
##
## Notes:
##
## Before running this script, run the data generating scripts:
##  * generate_linear_stochastic_data.R
##  * generate_driven_coupled_logistic_maps_data.R
##  * generate_lorenz_harmonic_oscillator_data.R
##  * generate_lorenz96_harmonic_oscillator_data.R
##
## ---------------------------

library(readr)
library(ggplot2)
library(patchwork)

source("R/plot_utils.R")
source("R/bootstrap.R")

#
# Linear stochastic system
#

lin_extreme <- read_csv(
  "data/linear_stochastic_extreme_coupling.csv",
  show_col_types = FALSE
)

plot_lin_extreme <- sync_plot(lin_extreme) +
  # scale_x_continuous(breaks = c(0, 0.02, 0.04, 0.06, 0.08)) +
  xlab(expression(coupling ~ strength ~ k))

#
# Coupled logistic maps with external driver
#

log_extreme <- read_csv(
  "data/driven_coupled_logistic_maps_extreme_coupling.csv",
  show_col_types = FALSE
) |>
  filter(RR_target %in% c(6, 10))

plot_log_extreme <- sync_plot(log_extreme) +
  # scale_x_continuous(breaks = c(0, 0.02, 0.04, 0.06, 0.08)) +
  xlab(expression(coupling ~ strength ~ eta))

#
# Lorenz-harmonic oscillator
#

lor_extreme <- read_csv(
  "data/lorenz_harmonic_extreme_coupling.csv",
  show_col_types = FALSE
)

plot_lor_extreme <- sync_plot(lor_extreme) +
  xlab(expression(coupling ~ strength ~ c))

#
# Lorenz 96  & harmonic oscillator
#

lor96_extreme <- read_csv(
  "data/lorenz96_harmonic_extreme_coupling.csv",
  show_col_types = FALSE
)

plot_lor96_extreme <- sync_plot(lor96_extreme) +
  xlab(expression(coupling ~ strength ~ kappa))

#
# Combined plot
#

extreme_plot <-
  plot_lin_extreme +
  plot_log_extreme + ylab("") +
  plot_lor_extreme + ylab("") +
  plot_lor96_extreme + ylab("") +
  plot_layout(nrow = 1) +
  plot_annotation(tag_levels = "A")

ggsave(extreme_plot, filename = "Plots/extreme_coupling.pdf",
       width = 25, height = 6)
