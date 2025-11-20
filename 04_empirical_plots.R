## ---------------------------
##
## Script name: 04_empirical_plots.R
##
## Purpose of script: Plot results of MvJRQA analysis and example time series
##
## Author: Dan Moenster & Sebastian Wallot
##
## Date Created: 6 November 2025
##
## ---------------------------
##
## Notes:
##
## The random seed is set to ensure reproducibility of bootstrapped
## error bars in the plots. Any deviation from the published results
## with a different random seed should be small enough to be ignored.
##
## ---------------------------

# Wipe work space
rm(list = ls())

library(readr)
library(edf)
library(tidyr)
library(ggplot2)
library(patchwork)
source("R/plot_utils.R")

set.seed(52328)

results <- read_csv("results/mvjrqa.csv", show_col_types = FALSE)

rr_long <- results |> 
  pivot_longer(
    cols = c(eye2d_JRR, eye3d_JRR,
             eye2d_DET, eye3d_DET,
             eye2d_LAM, eye3d_LAM,
             eye2d_DENTR, eye3d_DENTR,
             eye2d_RR1, eye3d_RR1,
             eye2d_RR2, eye3d_RR2,
             eye2d_RAD1, eye3d_RAD1,
             eye2d_RAD2, eye3d_RAD2,
             eye2d_file, eye3d_file),
    names_to = c("Eye_Type", ".value"),
    names_pattern = "eye(\\d+d)_(.*)"
  )


dimension_colours <- c(
  "2" = "#fd8d3c",
  "3" = "#b10026"
)

# Remove rows where RR1 deviates from target_RR by more than tolerance
jrci_plots <- list()
plot_count <- 1
tolerances <- c(0.2, 0.5, 1, 2, 4, 9)

for (tol in tolerances) {
  abs_rr_tol <- tol
  rr_jrci <- rr_long |>
    select(target_RR, Eye_Type, RR1, RR2, JRR) |> 
    filter(abs(RR1 - target_RR) <= abs_rr_tol) |>  # Abs. deviation
    rename(RR_target = target_RR) |> 
    mutate(coupling = if_else(Eye_Type == "2d", 2, 3)) |> 
    pivot_longer(cols = c("RR1", "RR2", "JRR"),
                 values_to = "RR",
                 names_to = "RR_type") |> 
    mutate(system = case_match(
      RR_type,
      "RR1" ~ "1",
      "RR2" ~ "2",
      "JRR" ~ "j"))
  
  ylimits <- c(0.008, 0.011)
  
  eeg_jrc <- jrc_plot(rr_jrci,
                      errorbars = "CI",
                      ci_level = 0.95,
                      samples  = 1e3,
                      coupling_name = "d") +
    scale_color_manual(values = dimension_colours,
                       aesthetics = c("colour", "fill")) +
    annotate("text", x = 6, y = .0110,
             label =  paste("delta == ", tol),
             parse = TRUE,
             size = 22 / .pt) +
    scale_y_continuous(labels = function(x) x * 100,
                       name = bquote(JRR/RR^2%*%100~"%"),
                       limits = ylimits)  +
    theme(text = element_text(size = 22))
  
  
  
  jrci_plots[[plot_count]] <- eeg_jrc
  plot_count <- plot_count + 1
}

eeg_jrc_delta <- patchwork::wrap_plots(jrci_plots) +
  plot_layout(axis_titles = "collect")

plot(eeg_jrc_delta)

ggsave(eeg_jrc_delta,
       filename = "plots/empirical_jrci_delta.pdf",
       width = 3 * 4, height = 3 * 3)

# Plot histograms of RR, indicating excluded observations.
hist_data <- data.frame()

for (tol in tolerances) {
  abs_rr_tol <- tol
  hist_data <- bind_rows(
    hist_data,
    rr_long |> 
      filter(Eye_Type == "2d") |> 
      mutate(delta = abs_rr_tol) |> 
      mutate(delta_str = gsub("Inf", "infinity", 
                              as.character(abs_rr_tol))) |> 
      mutate(excluded = if_else(abs(RR1 - target_RR) <= abs_rr_tol,
                                FALSE, TRUE)) |> 
      mutate(excluded_str = if_else(excluded,
                                    "excluded", "included"))
  )
}

hist_excluded <- hist_data |> 
  group_by(target_RR, delta) |> 
  summarise(n_removed = sum(excluded),
            n_total = n())

exclude_colors <- c(
  "excluded" = "#d95f02",
  "included" = "#1b9e77"
)

# Make custom labeller functions for the facet_grid
# Based on ggplot2's label_parsed
delta_labeller <- function (labels, multi_line = TRUE) 
{
  labels_raw <- label_value(labels, multi_line = multi_line)
  labels <- lapply(labels_raw, function(values) {
    paste0("delta == ", values)
  })
  if (multi_line) {
    lapply(unname(labels), lapply, function(values) {
      c(parse(text = as.character(values)))
    })
  }
  else {
    lapply(labels, function(values) {
      values <- paste0("list(", values, ")")
      lapply(values, function(expr) c(parse(text = expr)))
    })
  }
}

class(delta_labeller) <- c("function", "labeller")

rr_labeller <- function (labels, multi_line = TRUE) 
{
  labels_raw <- label_value(labels, multi_line = multi_line)
  labels <- lapply(labels_raw, function(values) {
    paste0("RR[T] == ", values)
  })
  if (multi_line) {
    lapply(unname(labels), lapply, function(values) {
      c(parse(text = as.character(values)))
    })
  }
  else {
    lapply(labels, function(values) {
      values <- paste0("list(", values, ")")
      lapply(values, function(expr) c(parse(text = expr)))
    })
  }
}

class(rr_labeller) <- c("function", "labeller")

delta_hist <- ggplot(hist_data,
                     aes(x = RR1 + 0.1, fill = excluded_str)) +
  geom_col(data = hist_excluded,
           aes(x = target_RR, y = n_total),
           width = 0.1,
           fill = NA,
           colour = "darkgray",
           inherit.aes = FALSE) +
  geom_histogram(alpha = 1, position = "identity", bins = 50, drop = TRUE) +
  geom_text(data = hist_excluded,
            aes(x = 2, y = 200, label = paste0("N[E] == ", n_removed)),
            parse = TRUE,
            inherit.aes = FALSE,
            size = 16 / .pt) +
  scale_fill_manual(values = exclude_colors,
                    guide = NULL) +
  scale_x_log10(name = bquote(RR),
                breaks = c(1, 4, 16, 64)) +
  scale_y_log10() +
  labs(y = "Frequency") +
  theme_classic() +
  theme(text = element_text(size = 22)) +
  facet_grid(target_RR ~ delta,
             labeller = labeller(
               delta= delta_labeller,
               target_RR = rr_labeller))

# plot(delta_hist)

ggsave(delta_hist,
       filename = "plots/empirical_delta_hist.pdf",
       width = 3 * 4, height = 3 * 4)


# JRCI plot for  delta = 1.5
abs_rr_tol <- 1.5
rr_jrci <- rr_long |>
  select(target_RR, Eye_Type, RR1, RR2, JRR) |> 
  filter(abs(RR1 - target_RR) <= abs_rr_tol) |> 
  rename(RR_target = target_RR) |> 
  mutate(coupling = if_else(Eye_Type == "2d", 2, 3)) |> 
  pivot_longer(cols = c("RR1", "RR2", "JRR"),
               values_to = "RR",
               names_to = "RR_type") |> 
  mutate(system = case_match(
    RR_type,
    "RR1" ~ "1",
    "RR2" ~ "2",
    "JRR" ~ "j"))

single_jrc <- jrc_plot(rr_jrci,
                    errorbars = "CI",
                    ci_level = 0.95,
                    samples  = 1e3,
                    coupling_name = "d") +
  scale_color_manual(values = dimension_colours,
                     aesthetics = c("colour", "fill")) +
  scale_y_continuous(labels = function(x) x * 100,
                     name = bquote(JRR/RR^2%*%100~"%"),
                     limits = c(0.93, 1.08) / 100)  + 
  theme(text = element_text(size = 22))


# plot(single_jrc)
ggsave(single_jrc,
       filename = "plots/empirical_jrci.pdf",
       width = 1.5 * 4, height = 1.5 * 3)

# Plot time series for a single trial


eeg <- read.edf("empirical_data/EEG_FLS/1_1_2.edf")
eeg <- data.frame(eeg$signal)
# Save the time from the second column, down sampled
eeg_time <- eeg[seq(1, dim(eeg)[1], by = 10), 2]
# remove t-variables and down sample to 50 Hz
eeg <- eeg[seq(1, dim(eeg)[1], by = 10), seq(1, dim(eeg)[2], by = 2)]
drops <- c("EEGHEOGRCPz.data",
           "EEGHEOGLCPz.data",
           "EEGVEOGUCPz.data",
           "EEGVEOGLCPz.data") # list EOGs to drop from data frame
eeg <- eeg[, !(names(eeg) %in% drops)] 
# z-score
eeg <- scale(eeg)
# Add the time as the first column
eeg <- cbind(eeg_time, eeg)
n_col <- ncol(eeg)
channels <- paste0("c", 1:(n_col - 1))
c_names <- c("t", channels)
colnames(eeg) <- c_names


# Eye tracking
df_eye <- read.csv("empirical_data/EYE_FLS/1_1_2.csv")
eye_2d <- data.frame(df_eye[, 1:2]) # extract 2d eye gaze
eye_3d <- data.frame(df_eye[, 3:5]) # extract 3d gaze
# cut off superfluous data points
if (dim(eeg)[1] > dim(df_eye)[1]) {
  eeg <- eeg[1:dim(eye_2d)[1], ] # cut eeg record
} else if (dim(eeg)[1] < dim(df_eye)[1]) {
  eye_2d <- eye_2d[1:dim(eeg)[1], ]# cut eye records
  eye_3d <- eye_3d[1:dim(eeg)[1], ]# cut eye records
} 

# scale() turned eeg into a matrix, now we need a data frame
eeg <- as.data.frame(eeg)
# Add time to eye tracking
eye_2d <- cbind(eeg$t, eye_2d)
colnames(eye_2d) <- c("t", "x", "y")

# Convert EEG to long format for easier plotting
eeg_long <- eeg |> 
  pivot_longer(cols = starts_with("c"),
               names_to = "channel",
               values_to = "potential") |> 
  mutate(channel = factor(channel,
                          levels = channels))

driver_colour <- "#1b9e77"
driven_colour <- "#7570b3"

eeg_plot <- ggplot(eeg_long,
                   aes(x = t, y = potential)) +
  geom_line(colour = driver_colour) +
  labs(
    x = "t",
    y = "Potential (normalized)"
  ) +
  scale_x_continuous(
    breaks = c(0, 100),
    guide = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(
    breaks = c(-2, 2),
    guide = guide_axis(minor.ticks = TRUE)) +
  theme_classic() +
  facet_wrap(~ channel, ncol = 14) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    text = element_text(size = 12)
  )

# plot(eeg_plot)

# Plot the path of fixations
xy_plot <- ggplot(eye_2d, aes(x = x, y = y)) +
  geom_path(colour = driven_colour) +
  coord_fixed(ratio = 1) +
  theme_classic()


# Convert eyue tracking to long format for easier plotting
eye_2d_long <- eye_2d |> 
  pivot_longer(cols = c("x", "y"),
               names_to = "variable",
               values_to = "value") 

x_plot <- ggplot(eye_2d, aes(x = t, y = x)) +
  geom_line(colour = driven_colour, linewidth = 1) +
  labs(
    x = "t",
    y = "x"
  ) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

y_plot <- ggplot(eye_2d, aes(x = t, y = y)) +
  geom_line(colour = driven_colour, linewidth = 1) +
  labs(
    x = "t",
    y = "y"
  ) +
  theme_classic()

eye_plot <- (x_plot / y_plot | xy_plot) &
  theme(text = element_text(size = 12),
        axis.title.y = element_text(angle = 0,  vjust = 0.5, hjust = 1))

# plot(eye_plot)

all_plot <- eeg_plot / eye_plot +
  plot_layout(heights = c(2.5, 1)) +
  plot_annotation(tag_levels = list(c("A", "B", "", "C")))

# plot(all_plot)

ggsave("plots/empirical_time_series_plot.pdf",
       plot = all_plot,
       width = 7.5, height = 10)

