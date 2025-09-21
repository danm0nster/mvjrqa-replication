## ---------------------------
##
## Script name: plot_high_dimensional_lorenz96.R
##
## Purpose of script: Plot results high-dimensional Lorenz 96 system.
##
## Author: Dan Moenster & Sebastian Wallot
##
## Date Created: Jun 4 2025
##
## ---------------------------
##
## Notes:
##
## This script produces plots for the high-dimensional Lorenz 96 system
## coupled to a harmonic oscillator.
##
## It uses data produced by the data generation script
## generate_lorenz96_high_dimensional_oscillator_data.R.
##
## ---------------------------

source("R/lorenz96.R")
source("R/plot_utils.R")
source("R/utils.R")

# Make sure "./plots" directory exists
create_dir_if_not_present("plots")

#
# Time series plots
#
n_var <- 16
model_data <- lorenz96_harmonic(n = 1000,
                                skip = 1000,
                                n_var = n_var,
                                coupling = 0.5)

plot_data <- model_data |>
  pivot_longer(cols = -time,
               values_to = "value",
               names_to = "var") |>
  # Create expression to use with plotmath and the label_parsed labeller
  mutate(var_expression = gsub("x([0-9]+)", "x[\\1]", var)) |>
  mutate(var_expression = factor(var_expression,
                                 levels = c(paste0("x[", seq(1, n_var), "]"),
                                            "u", "v"))) |>
  mutate(system = if_else(var %in% c("u", "v"), 2, 1)) |>
  mutate(system = factor(system))

system_colours <- c(
  "1" = "#1b9e77",
  "2" = "#7570b3"
)

# Only include time label and ticks on the x-axis on the last (bottom) plot
lor96_hd_plot <- ggplot(plot_data, aes(x = time, y = value, colour = system)) +
  geom_line(linewidth = 2) +
  xlab("t") +
  ylab(expression("")) +
  scale_x_continuous(breaks = seq(10, 20, by = 2)) +
  scale_y_continuous(breaks = scales::breaks_pretty(3)) +
  scale_colour_manual(values = system_colours,
                      aesthetics = c("colour", "fill")) +
  guides(colour = "none", fill = "none") +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0,  vjust = 0.5, hjust = 1)
  ) +
  facet_wrap(~ var_expression, ncol = 2,
             strip.position = "left",
             axes = "all_x",
             axis.labels = "margins",
             scales = "free_y",
             labeller = label_parsed) +
  theme(text = element_text(size = 28))

ggsave("plots/lorenz96_high_dim_time_series_plot.pdf",
       plot = lor96_plot,
       width = 10, height = 20)
#
# High-dimensional Lorenz 96 system coupled to harmonic oscillator
#

lorenz96_high_dim_coupling_sweep <- read_csv(
  "data/lorenz96_high_dim_coupling_sweep.csv"
)

lor96_hid_jrr_rr_plot <- sync_plot(lorenz96_high_dim_coupling_sweep) +
  xlab(expression(coupling ~ strength ~ kappa)) +
  expand_limits(x = 0.505)

ggsave(plot = lor96_hid_jrr_rr_plot,
       filename = "plots/lorenz96_high_dim_harmonic_jrr_rr.pdf",
       width = 1.5 * 4, height = 1.5 * 3)


lorenz_coupling_colours <- c(
  "0" = "#fd8d3c",
  "0.1" = "#fc4e2a",
  "0.2" = "#e31a1c",
  "0.4" = "#b10026"
)

lorenz_hi_dimc_rr_sweep <- read_csv(
  "data/lorenz96_high_dim_rr_sweep.csv"
)



lor_hi_dim_jrc_plot <- jrc_plot(lorenz_hi_dimc_rr_sweep,
                                errorbars = "CI",
                                ci_level = 0.95,
                                coupling_name = "kappa") +
  scale_color_manual(values = lorenz_coupling_colours,
                     aesthetics = c("colour", "fill")) +
  expand_limits(x = 0, y = 1.4) +
  scale_x_log10(breaks = c(0.1, 1, 10, 100), labels = c(0.1, 1, 10, 100)) +
  scale_y_continuous(breaks = seq(0, 1.25, by = 0.25),
                     labels = c("0.00", "0.25", "0.50",
                                "0.75", "1.00", "1.25")) +
  theme(text = element_text(size = 22))

ggsave(lor_hi_dim_jrc_plot,
       filename = "plots/lor_hi_dim_jrr_by_rr2_log_color.pdf",
       width = 1.5 * 4, height = 1.5 * 3)


#
# Extreme coupling
#

lor96_hi_dim_extreme <- read_csv("data/lorenz96_high_dim_extreme_coupling.csv")

plot_lor96_hi_dim_extreme <- sync_plot(lor96_hi_dim_extreme) +
  xlab(expression(coupling~strength~kappa)) +
  expand_limits(x = 1.505) # The last x-label at 1.5 gets cropped, so add space

plot(plot_lor96_hi_dim_extreme)

ggsave(plot = plot_lor96_hi_dim_extreme,
       filename = "plots/lorenz96_high_dim_extreme_coupling.pdf",
       width = 1.5 * 4, height = 1.5 * 3)
