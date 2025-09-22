## ---------------------------
##
## Script name: plot_model_time_series.R
##
## Purpose of script: Plot example time series for the model systems.
##
## Author: Dan Moenster & Sebastian Wallot
##
## Date Created: Jun 17 2025
##
## ---------------------------
##
## Notes:
## The plots use relatively short time series, so the data are not saved,
## but rather recomputed when the script is run.
##
## ---------------------------

library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)

source("R/coupled_logistic_map.R")
source("R/lorenz_harmonic.R")
source("R/lorenz96.R")
source("R/utils.R")

driver_colour <- "#1b9e77"
driven_colour <- "#7570b3"

#
# Set random seed for reproducibility
#
set.seed(85322)

#
# Plot example of the stochastic linear model
#
b_1 <- 2
sample_size <- 100
epsilon_4 <- rnorm(sample_size)
x_1 <- rnorm(sample_size)
y_1 <- b_1 * x_1 + epsilon_4 + rnorm(sample_size)
y_2 <- b_1 * x_1 + epsilon_4 + rnorm(sample_size)

lin_data <- data.frame(
  x_1 = x_1,
  y_1 = y_1,
  y_2 = y_2,
  time = 1:sample_size
)

x_1_plot <- ggplot(lin_data, aes(x = time, y = x_1)) +
  geom_line(colour = driver_colour, linewidth = 2) +
  ylab(expression(x)) +
  scale_y_continuous(breaks = c(-2, 0, 2)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

y_1_plot <- ggplot(lin_data, aes(x = time, y = y_1)) +
  geom_line(colour = driven_colour, linewidth = 2) +
  ylab(expression(y[1])) +
  scale_y_continuous(breaks = c(-4, 0, 4)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

y_2_plot <- ggplot(lin_data, aes(x = time, y = y_2)) +
  geom_line(colour = driven_colour, linewidth = 2) +
  xlab("t") +
  ylab(expression(y[2])) +
  scale_y_continuous(breaks = c(-4, 0, 4)) +
  theme_classic()

all_plot <- x_1_plot / y_1_plot / y_2_plot &
  theme(text = element_text(size = 22),
        axis.title.y = element_text(angle = 0,  vjust = 0.5, hjust = 1))

ggsave("plots/linear_model_time_series_plot.pdf",
       plot = all_plot,
       width = 6, height = 7)


#
# Plot of the coupled logistic map with external driver
#

x_init <- runif(1, min = 0.2, max = 0.8)
y_init <- runif(1, min = 0.2, max = 0.8)
phi_init <- runif(1, min = 0, max = 2 * pi)
rx <-  3.58
ry <-  3.6
byx <- 0.1
bxy <-  0
cc <- 0.2

map_data <- coupled_logistic_map(
  x0 = x_init,
  y0 = y_init,
  rx = rx,
  ry = ry,
  bxy = bxy,
  byx = byx,
  etax = cc,
  etay = cc,
  p = 30,
  phi = phi_init,
  N = 100,
  N_skip = 1000
)

map_h_plot <- ggplot(map_data, aes(x = time, y = H)) +
  geom_line(colour = driver_colour, linewidth = 2) +
  xlab("t") +
  ylab("H") +
  scale_y_continuous(breaks = c(-1, 0, 1)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

map_x_plot <- ggplot(map_data, aes(x = time, y = X)) +
  geom_line(colour = driven_colour, linewidth = 2) +
  xlab("t") +
  ylab("x") +
  scale_y_continuous(breaks = c(0.3, 0.6, 0.9)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

map_y_plot <- ggplot(map_data, aes(x = time, y = Y)) +
  geom_line(colour = driven_colour, linewidth = 2) +
  xlab("t") +
  ylab("y") +
  scale_y_continuous(breaks = c(0.3, 0.6, 0.9)) +
  theme_classic()


map_plot <- map_h_plot / map_x_plot / map_y_plot &
  theme(text = element_text(size = 22),
        axis.title.y = element_text(angle = 0,  vjust = 0.5, hjust = 1))

ggsave("plots/logistic_driver_time_series_plot.pdf",
       plot = map_plot,
       width = 6, height = 7)


#
# Example time series for the Lorenz system and harmonic oscillator
#

data_length <- 1000

model_data <- lorenz_harmonic(n = data_length, skip = 1000, coupling = 0.4)

# Only include time label and ticks on the x-axis on the last (bottom) plot
lor_x_plot <- ggplot(model_data, aes(x = time, y = x)) +
  geom_line(colour = driver_colour, linewidth = 2) +
  xlab("t") +
  ylab("x") +
  scale_y_continuous(breaks = c(-15, 0, 15)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

lor_y_plot <- ggplot(model_data, aes(x = time, y = y)) +
  geom_line(colour = driver_colour, linewidth = 2) +
  xlab("t") +
  ylab("y") +
  scale_y_continuous(breaks = c(-15, 0, 15)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

lor_z_plot <- ggplot(model_data, aes(x = time, y = z)) +
  geom_line(colour = driver_colour, linewidth = 2) +
  xlab("t") +
  ylab("z") +
  scale_y_continuous(breaks = c(5, 25, 45)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

osc_u_plot <- ggplot(model_data, aes(x = time, y = u)) +
  geom_line(colour = driven_colour, linewidth = 2) +
  xlab("t") +
  ylab("u") +
  scale_y_continuous(breaks = c(-10, 0, 10)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  expand_limits(y = -10)

osc_v_plot <- ggplot(model_data, aes(x = time, y = v)) +
  geom_line(colour = driven_colour, linewidth = 2) +
  xlab("t") +
  ylab("v") +
  scale_x_continuous(breaks = seq(10, 20, by = 2)) +
  scale_y_continuous(breaks = c(-25, 0, 25)) +
  theme_classic()

lor_osc_plot <- lor_x_plot / lor_y_plot / lor_z_plot / osc_u_plot / osc_v_plot &
  theme(text = element_text(size = 22),
        axis.title.y = element_text(angle = 0,  vjust = 0.5, hjust = 1))

ggsave("plots/lorenz_harmonic_time_series_plot.pdf",
       plot = lor_osc_plot,
       width = 6, height = 7)

#
# Lorenz 96 & harmonic oscillator time series plots
#

l96_data <- lorenz96_harmonic(n = 1000, skip = 1000, coupling = 0.8)

l96_plot_data <- l96_data |>
  pivot_longer(cols = -time,
               values_to = "value",
               names_to = "var") |>
  # Create expression to use with plotmath and the label_parsed labeller
  mutate(var_expression = gsub("x([1-5]+)", "x[\\1]", var)) |>
  mutate(var_expression = factor(var_expression,
                                 levels = c(paste0("x[", seq(1, 5), "]"),
                                            "u", "v"))) |>
  mutate(system = if_else(var %in% c("u", "v"), 2, 1)) |>
  mutate(system = factor(system))

system_colours <- c(
  "1" = driver_colour,
  "2" = driven_colour
)

custom_breaks <- function(x) {
  # We want symmetric scales around 0, although the data might be skewed.
  # We do not have room for many breaks, so keep it simple.
  smallest_endpoint <- floor(min(abs(min(x)), abs(max(x))))
  # We want either (-5, 5), (-20, 20) or (-40, 40)
  breaks_vector <- c(5, 15, 40)
  a <- breaks_vector[which.min(abs(breaks_vector - smallest_endpoint))]
  c(-a, 0, a)
}

# Only include time label and ticks on the x-axis on the last (bottom) plot
lor96_plot <- ggplot(l96_plot_data, aes(x = time, y = value, colour = system)) +
  geom_line(linewidth = 2) +
  xlab("t") +
  ylab(expression("")) +
  scale_x_continuous(breaks = seq(10, 20, by = 2)) +
  scale_y_continuous(breaks = custom_breaks) +
  scale_colour_manual(values = system_colours,
                      aesthetics = c("colour", "fill")) +
  guides(colour = "none", fill = "none") +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0,  vjust = 0.5, hjust = 1)) +
  facet_wrap(~ var_expression, nrow = 7,
             strip.position = "left",
             axes = "all_x",
             axis.labels = "margins",
             scales = "free_y",
             labeller = label_parsed) +
  theme(text = element_text(size = 22))

ggsave("plots/lorenz96_time_series_plot.pdf",
       plot = lor96_plot,
       width = 6, height = 7)
