## ---------------------------
##
## Script name: plot_phase_space_coupling.R
##
## Purpose of script: Plot coupling between Lorenz and harmonic oscillator
##
## Author: Dan Moenster & Sebastian Wallot
##
## Date Created: Aug 22 2023
##
## ---------------------------
##
## Notes:
##
## This script creates a plot that shows how the attractor of the harmonic
## oscillator changes with coupling strength from the Lorenz system.
##
## ---------------------------

library(ggplot2)
library(patchwork)

source("R/lorenz_harmonic.R")

# None of the time series need to be further embedded,
# so delay and embed are set to 1.
delay <- 1
embed <- 1

# Set radius for time series 1
radius_1 <- 0.2
# Radius for time series 2 should be fixed to give same RR, but here we
# just do a conceptual example, so we set it to the same value as ts 1.
radius_2 <- 0.2

data_length <- 2000

plot_list <- list()
c_values <- c(0, 0.1, 0.2)
plot_num <- 0
for (c in c_values) {
  plot_num <- plot_num + 1
  model_data <- lorenz_harmonic(n = data_length, skip = 1000, coupling = c)
  plot_list[[plot_num]] <- ggplot(model_data, aes(x = u, y = v)) +
    geom_path(alpha = 0.5) + 
    theme_classic() +
    ggtitle(paste("c = ", c)) +
    theme(text = element_text(size = 14))
}

for (k in 1:length(plot_list)) {
  if (k == 1) {
    figure <- plot_list[[k]]
  } else {
    figure <- figure + plot_list[[k]]
  }
}

# print(figure)

ggsave(figure, filename = "Plots/harmonic_phase_space_coupling.pdf", 
       width = 9, height = 3)


