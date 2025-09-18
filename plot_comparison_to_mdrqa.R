## ---------------------------
##
## Script name: plot_comparison_to_mdrqa.R
##
## Purpose of script: Produce plots comparing MvJRQA to MdRQA
##
## Author: Dan Moenster & Sebastian Wallot
##
## Date Created: May 18 2025
##
## ---------------------------
##
## Notes:
##
## Before running this script, run the data generating script:
## generate_data_for_comparison_to_mdrqa.R
##
## ---------------------------

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

source("R/utils.R")

#
# Linear stochastic system
#

linear_mdrqa_comparison <- read_csv(
  "data/linear_stochastic_system_mdrqa_comparison.csv"
)

# Convert to long format first
linear_plot_data <- linear_mdrqa_comparison |>
  select(index_simu, coupling, jrqa_RR, mdrqa_RR) |>
  rename(sample_id = index_simu) |>
  pivot_longer(cols = ends_with("_RR"),
               values_to = "RR",
               names_to = "system") |>
  mutate(system = if_else(system == "jrqa_RR", "j", "f"))

linear_summary_data <- linear_plot_data |>
  filter(system %in% c("j", "f")) |>
  group_by(coupling, system) |>
  summarise(RR_M = mean(RR),
            RR_SE = sd(RR) / sqrt(n()),
            CI_lo = lower_ci(RR),
            CI_hi = upper_ci(RR)) |>
  ungroup()


#
# Coupled logistic map with external driver
#

logistic_mdrqa_comparison <- read_csv(
  "data/logistic_system_mdrqa_comparison.csv"
)

logistic_plot_data <- left_join(
  logistic_mdrqa_comparison |>
    filter(system == "j") |>
    select(sample_id, coupling, RR) |>
    rename(RR_j = RR),
  logistic_mdrqa_comparison |>
    filter(system == "f") |>
    select(sample_id, coupling, RR) |>
    rename(RR_m = RR),
)

logistic_summary_data <- logistic_mdrqa_comparison |>
  filter(system %in% c("j", "f")) |>
  group_by(coupling, system) |>
  summarise(RR_M = mean(RR),
            RR_SE = sd(RR) / sqrt(n()),
            CI_lo = lower_ci(RR),
            CI_hi = upper_ci(RR)) |>
  ungroup()

#
# Lorenz and harmonic oscillator
#

lorenz_harmonic_mdrqa_comparison <- read_csv(
  "data/lorenz_harmonic_mdrqa_comparison.csv"
)

lorenz_plot_data <- left_join(
  lorenz_harmonic_mdrqa_comparison |>
    filter(system == "j") |>
    select(sample_id, coupling, RR) |>
    rename(RR_j = RR),
  lorenz_harmonic_mdrqa_comparison |>
    filter(system == "f") |>
    select(sample_id, coupling, RR) |>
    rename(RR_m = RR),
)

lorenz_summary_data <- lorenz_harmonic_mdrqa_comparison |>
  filter(system %in% c("j", "f")) |>
  group_by(coupling, system) |>
  summarise(RR_M = mean(RR),
            RR_SE = sd(RR) / sqrt(n()),
            CI_lo = lower_ci(RR),
            CI_hi = upper_ci(RR)) |>
  ungroup()

#
# Lorenz 96 and harmonic oscillator
#

lorenz96_harmonic_mdrqa_comparison <- read_csv(
  "data/lorenz_96_harmonic_mdrqa_comparison.csv"
)

lorenz96_plot_data <- left_join(
  lorenz96_harmonic_mdrqa_comparison |>
    filter(system == "j") |>
    select(sample_id, coupling, RR) |>
    rename(RR_j = RR),
  lorenz96_harmonic_mdrqa_comparison |>
    filter(system == "f") |>
    select(sample_id, coupling, RR) |>
    rename(RR_m = RR),
)

lorenz96_summary_data <- lorenz96_harmonic_mdrqa_comparison |>
  filter(system %in% c("j", "f")) |>
  group_by(coupling, system) |>
  summarise(RR_M = mean(RR),
            RR_SE = sd(RR) / sqrt(n()),
            CI_lo = lower_ci(RR),
            CI_hi = upper_ci(RR)) |>
  ungroup()

#
# Plot RR vs. coupling for MvJRQA and MdRQA
#

# Linear stochastic system

linear_label_data <- linear_summary_data |>
  filter(coupling == max(coupling)) |>
  mutate(label = if_else(system == "j", "MvJRQA", "MdRQA")) |>
  mutate(RR_M = if_else(system == "j", RR_M - 2, RR_M + 1))

errorbar_width <- 0.3 * max(linear_summary_data$coupling) /
  length(unique(linear_summary_data$coupling))

p_linear <- ggplot(linear_summary_data,
                   aes(x = coupling,
                       y = RR_M,
                       group = system,
                       colour = system)) +
  geom_errorbar(aes(ymin = CI_lo,
                    ymax = CI_hi),
                width = errorbar_width) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("j" = "#1f78b4", "f" = "#33a02c"),
                      guide = NULL) +
  geom_text(data = linear_label_data,
            aes(label = label),
            parse = TRUE,
            nudge_x = 0,
            hjust = "right",
            size = 5) +
  xlab(expression(Coupling ~ k)) +
  ylab("RR") +
  theme_classic() +
  theme(text = element_text(size = 22))

plot(p_linear)


# Logistic maps

logistic_label_data <- logistic_summary_data |>
  filter(coupling == min(coupling)) |>
  mutate(label = if_else(system == "j", "MvJRQA", "MdRQA")) |>
  mutate(RR_M = if_else(system == "j", RR_M - 0.08, RR_M + 0.08))

errorbar_width <- 0.3 * max(logistic_summary_data$coupling) /
  length(unique(logistic_summary_data$coupling))

p_logistic <- ggplot(logistic_summary_data,
                     aes(x = coupling,
                         y = RR_M,
                         group = system,
                         colour = system)) +
  geom_errorbar(aes(ymin = CI_lo,
                    ymax = CI_hi),
                width = errorbar_width) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("j" = "#1f78b4", "f" = "#33a02c"),
                      guide = NULL) +
  geom_text(data = logistic_label_data,
            aes(label = label),
            parse = TRUE,
            nudge_x = 0,
            hjust = "left",
            size = 5) +
  xlab(expression(Coupling ~ eta)) +
  ylab("RR") +
  theme_classic() +
  theme(text = element_text(size = 22))

plot(p_logistic)

# Lorenz-harmonic oscillator

lorenz_label_data <- lorenz_summary_data |>
  filter(coupling == max(coupling)) |>
  mutate(label = if_else(system == "j", "MvJRQA", "MdRQA")) |>
  mutate(RR_M = if_else(system == "j", RR_M + 0.03, RR_M - 0.02))

errorbar_width <- 0.3 * max(lorenz_summary_data$coupling) /
  length(unique(lorenz_summary_data$coupling))

p_lorenz <- ggplot(lorenz_summary_data,
                   aes(x = coupling,
                       y = RR_M,
                       group = system,
                       colour = system)) +
  geom_errorbar(aes(ymin = CI_lo,
                    ymax = CI_hi),
                width = errorbar_width) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("j" = "#1f78b4", "f" = "#33a02c"),
                      guide = NULL) +
  geom_text(data = lorenz_label_data,
            aes(label = label),
            parse = TRUE,
            nudge_x = 0,
            hjust = "right",
            size = 5) +
  xlab(expression(Coupling ~ c)) +
  ylab("RR") +
  theme_classic() +
  theme(text = element_text(size = 22))

plot(p_lorenz)

# Lorenz 96  & harmonic oscillator

lorenz96_label_data <- lorenz96_summary_data |>
  filter(coupling == max(coupling)) |>
  mutate(label = if_else(system == "j", "MvJRQA", "MdRQA")) |>
  mutate(RR_M = if_else(system == "j", RR_M + 0.03, RR_M - 0.02))

errorbar_width <- 0.3 * max(lorenz96_summary_data$coupling) /
  length(unique(lorenz96_summary_data$coupling))

p_lorenz96 <- ggplot(lorenz96_summary_data,
                     aes(x = coupling,
                         y = RR_M,
                         group = system,
                         colour = system)) +
  geom_errorbar(aes(ymin = CI_lo,
                    ymax = CI_hi),
                width = errorbar_width) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("j" = "#1f78b4", "f" = "#33a02c"),
                      guide = NULL) +
  geom_text(data = lorenz96_label_data,
            aes(label = label),
            parse = TRUE,
            nudge_x = 0,
            hjust = "right",
            size = 5) +
  xlab(expression(Coupling ~ kappa)) +
  ylab("RR") +
  theme_classic() +
  theme(text = element_text(size = 22))

plot(p_lorenz96)

comparison_plot <-  p_linear + p_logistic + p_lorenz + p_lorenz96 +
  plot_annotation(tag_levels = "A") +
  plot_layout(nrow = 1)

plot(comparison_plot)
create_dir_if_not_present("plots")
ggsave("plots/comparison_MdRQA_MvJRQA_coupling.pdf", width = 20, height = 5)
