#
# Calculate data for the Lorenz system coupled to a harmonic oscillator
#

source("R/data_sampling.R")

# Sweep over moderate coupling values
rr_values <- c(1, 2, 5, 10)
coupling_values <- seq(0, 0.3, 0.025)

# Set random seed to make the script reproducible
set.seed(64225)

compute_lorenz_harmonic_data(
  coupling = coupling_values,
  rr = rr_values,
  num_samples = 100,
  sample_length = 500,
  sample_skip = 100,
  csv_file_name = "lorenz_harmonic_coupling_sweep.csv"
)

# Extreme coupling
rr_values <- c(1, 2, 5, 10)
coupling_values_extreme <- seq(0, 1.2, 0.05)

# Set random seed to make the script reproducible
set.seed(87435)

compute_lorenz_harmonic_data(
  coupling = coupling_values_extreme,
  rr = rr_values,
  num_samples = 100,
  sample_length = 500,
  sample_skip = 100,
  csv_file_name = "lorenz_harmonic_extreme_coupling.csv"
)

# RR sweep
coupling_values <- c(0, 0.1, 0.2, 0.4)
rr_values <- c(0.2, 0.3, 0.6, 1.25, 2.5, 5, 10, 20, 40, 80)

# Set random seed to make the script reproducible
set.seed(64382)

compute_lorenz_harmonic_data(
  coupling = coupling_values,
  rr = rr_values,
  num_samples = 100,
  sample_length = 500,
  sample_skip = 100,
  csv_file_name = "lorenz_harmonic_rr_sweep.csv"
)

