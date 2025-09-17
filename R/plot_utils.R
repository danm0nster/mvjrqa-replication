## ---------------------------
##
## File name: plot_utils.R
##
## Purpose:  Utility functions for plotting
##
## Author: Dan Moenster & Sebastian Wallot
##
## Date Created: May 17 2025
##
## ---------------------------
##
## Notes:
##
## Function for plotting the Joint Recurrence Coupling Indicator: jrc_plot()
## Function for plotting synchrony (JRR / RR): sync_plot()
##
## ---------------------------

library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)
library(latex2exp)
library(ggrepel)


source("R/utils.R")
# TODO: The summary functions here could use the function
#       compute_sample_summary() from utils.R instead of replicating code here.

########### Function to produce Joint Recurrence Coupling Plot
jrc_plot <- function(joint_rr,
                     errorbars = c("CI", "SE", "Off"),
                     ci_level = 0.95, # Level for error bars when "CI"
                     coupling_name = "eta",
                     show_zero = FALSE, # If TRUE, will plot line at zero
                     plot_title = NULL,
                     csv_file = NULL) {
  errorbars <- match.arg(errorbars)
  joint_rr_summary <- joint_rr |>
    group_by(RR_target, coupling, system) |>
    summarise(
      rr_sd = sd(RR),
      rr_se = rr_sd / sqrt(n()),
      rr_mean = mean(RR),
      .groups = "drop"
    ) |>
    ungroup()
  # Find the largest ordinate value in the data to be plotted
  max_data_val <- joint_rr_summary |>
    filter(system == "j") |>
    mutate(jrr_div_rr2 = rr_mean / RR_target^2) |>
    pull(jrr_div_rr2) |>
    max()
  rr_identical_min <- 1 / max_data_val
  rr_random_min <- min(joint_rr$RR_target)
  N_groups <- nrow(joint_rr) %/% 3
  joint_rr_wide <- bind_cols(
    data.frame(group_id = rep(1:N_groups, each = 3)),
    joint_rr
  ) |>
    select(group_id, RR, system, coupling, RR_target) |>
    pivot_wider(names_from = "system", names_prefix = "RR_",
                values_from = "RR") |>
    mutate(RR_12_mean = 0.5 * (RR_1 + RR_2)) |>
    # Compute the Joint Recurrence Coupling Indicator (JRCI)
    # This is a tentative name for JRCI = JRR / (0.5 * (RR_1 + RR_2))^2
    #
    mutate(jrci = RR_j / RR_12_mean^2)
  # Make a summary with means, sd, se and CI
  joint_rr_summary_wide <- joint_rr_wide |>
    group_by(RR_target, coupling) |>
    summarise(
      rr_sd = sd(RR_j),
      rr_se = rr_sd / sqrt(n()),
      rr_mean = mean(RR_j),
      rr_sub_sd = sd(RR_12_mean),
      rr_sub_se = rr_sub_sd / sqrt(n()),
      rr_sub_mean = mean(RR_12_mean),
      jrci_mean = mean(jrci),
      jrci_sd = sd(jrci),
      jrci_se = jrci_sd / sqrt(n()),
      jrci_ci_lo = lower_ci(jrci),
      jrci_ci_hi = upper_ci(jrci),
      .groups = "drop"
    ) |>
    ungroup()
  # Build the plot
  coupling_plot <- ggplot(joint_rr_summary_wide,
                          aes(x = RR_target,
                              y = jrci_mean,
                              colour = factor(coupling),
                              fill = factor(coupling),
                              group = factor(coupling)))
  if (errorbars != "Off") {
    if (errorbars == "CI") {
      coupling_plot <- coupling_plot +
        geom_errorbar(aes(ymin = jrci_ci_lo,
                          ymax = jrci_ci_hi),
                      width = 0.05,
                      linetype = "solid")
    } else if (errorbars == "SE") {
      coupling_plot <- coupling_plot +
        geom_errorbar(aes(ymin = jrci_mean - jrci_se,
                          ymax = jrci_mean + jrci_se),
                      width = 0.05,
                      linetype = "solid")
    }
  }
  coupling_plot <- coupling_plot +
    geom_point(alpha = 0.8) +
    geom_line(linewidth = 2) +
    # Random null model
    geom_line(data = data.frame(RR_target = seq(rr_random_min, 100, 1)),
              aes(x = RR_target, y = 1 / 100),
              inherit.aes = FALSE,
              color = "blue") +
    # Identical systems model
    geom_line(data = data.frame(RR_target = seq(rr_identical_min, 100, 0.1)),
              aes(x = RR_target, y = 1 / RR_target),
              inherit.aes = FALSE,
              color = "blue") +
    scale_x_log10() +
    xlab("Sub-system RR (%)") +
    ylab(TeX("$JRR / RR^2$")) +
    # The nudge gives warnings because of infinities in the log-scale, but
    # seems to work, nonetheless.
    geom_text_repel(data = joint_rr_summary_wide |>
                      filter(RR_target == min(RR_target)),
                    aes(label = paste(coupling_name, " == ", coupling)),
                    parse = TRUE,
                    size = 5,
                    direction = "y", nudge_x = -0.5, hjust = "left",
                    segment.colour = NA) +
    theme_classic() +
    theme(legend.position = "none") +
    expand_limits(x = -0.5) + # To make room for geom_text_repel() labels
    expand_limits(x = 100)
  if (!is.null(plot_title)) {
    coupling_plot <- coupling_plot + labs(title = plot_title)
  }
  if (show_zero) {
    # Theoretical minimum is zero
    coupling_plot <- coupling_plot +
      geom_line(data = data.frame(RR_target = seq(rr_random_min, 100, 1)),
                aes(x = RR_target, y = 0),
                inherit.aes = FALSE,
                linetype = "dashed",
                color = "grey")
  }
  if (!is.null(csv_file)) {
    readr::write_csv(joint_rr_summary_wide, file = csv_file)
  }
  return(coupling_plot)
}

sync_plot <- function(joint_rr,
                      errorbars = c("CI", "SE", "Off"),
                      ci_level = 0.95,
                      csv_file = NULL) {
  errorbars <- match.arg(errorbars)
  min_coupling <- min(joint_rr$coupling)
  max_coupling <- max(joint_rr$coupling)
  error_bar_width <- (max_coupling - min_coupling) / 50
  # 1. Convert to wide format. Keep RR measure
  # 2. Compute synchrony = RR_j / 0.5 * (RR_1 + RR_2)
  # 3. Compute CIs for synchrony
  # 4. Make the plot
  # Convert to wide format and compute synchrony for each sample
  joint_rr_wide <-  joint_rr |>
    select(sample_id, RR, system, coupling, RR_target) |>
    pivot_wider(names_from = "system", names_prefix = "RR_",
                values_from = "RR") |>
    mutate(RR_12_mean = 0.5 * (RR_1 + RR_2)) |>
    mutate(synchrony = RR_j / RR_12_mean)
  # Make a summary with mean, SD, SD and CI
  joint_rr_summary_wide <- joint_rr_wide |>
    group_by(RR_target, coupling) |>
    summarise(
      rr_sd = sd(RR_j),
      rr_se = rr_sd / sqrt(n()),
      rr_mean = mean(RR_j),
      rr_sub_sd = sd(RR_12_mean),
      rr_sub_se = rr_sub_sd / sqrt(n()),
      rr_sub_mean = mean(RR_12_mean),
      sync_mean = mean(synchrony),
      sync_sd = sd(synchrony),
      sync_se = sync_sd / sqrt(n()),
      sync_ci_lo = lower_ci(synchrony),
      sync_ci_hi = upper_ci(synchrony),
      .groups = "drop"
    ) |>
    ungroup()
  s_plot <- ggplot(joint_rr_summary_wide,
                   aes(x = coupling, y = sync_mean,
                       group = RR_target,
                       colour = RR_target)) +
    scale_colour_gradient(low = "#56B1F7",
                          high = "#132B43",
                          aesthetics = c("colour", "fill"))
  if (errorbars != "Off") {
    if (errorbars == "CI") {
      s_plot <- s_plot +
        geom_errorbar(aes(ymin = sync_ci_lo,
                          ymax = sync_ci_hi),
                      width = error_bar_width,
                      linetype = "solid")
    } else if (errorbars == "SE") {
      s_plot <- s_plot +
        geom_errorbar(aes(ymin = sync_mean - sync_se,
                          ymax = sync_mean + sync_se),
                      width = error_bar_width,
                      linetype = "solid")
    }
  }
  s_plot <- s_plot +
    geom_point(size = 3) +
    geom_line(linewidth = 1.5) +
    geom_text_repel(data = joint_rr_summary_wide |>
                      filter(coupling == max(coupling)),
                    aes(label = paste("RR == ", RR_target, "*'%'")),
                    parse = TRUE,
                    direction = "y",
                    nudge_x = 0.2 * max(joint_rr$coupling),
                    hjust = "left",
                    size = 5,
                    segment.colour = NA) +
    xlab("coupling strength c") +
    ylab("JRR / RR") +
    theme_classic() +
    theme(legend.position = "none") +
    theme(text = element_text(size = 22)) +
    expand_limits(x = max(joint_rr$coupling) * 1.03)
  if (!is.null(csv_file)) {
    readr::write_csv(joint_rr_summary_wide, file = csv_file)
  }
  return(s_plot)
}
