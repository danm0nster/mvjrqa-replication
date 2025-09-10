# load libraries
library(crqa)
library(readr)

# Set random seed for replicability
set.seed(98879)

# source functions
source("R/mdrqa_comparison_sampling.R")


#
# Linear stochastic system
#

compute_linear_stochastic_mdrqa_comparison(
  coupling = seq(0, 10, by = 0.5),
  rr = 20,
  mdrqa_radius = 0.6,
  num_samples = 100,
  sample_length = 500,
  csv_file_name = "linear_stochastic_system_mdrqa_comparison.csv"
)

#
# Coupled logistic maps with external driver
#

compute_logistic_system_mdrqa_comparison(
  rr = 6,
  mdrqa_radius = 0.75,
  coupling = seq(0, 0.04, by = 0.04 / 20),
  num_samples = 100,
  sample_length = 500,
  sample_skip = 100,
  csv_file_name = "logistic_system_mdrqa_comparison.csv"
)

#
# Lorenz-Harmonic oscillator
#

compute_lorenz_harmonic_mdrqa_comparison(
  rr = 0.5,
  coupling = seq(0, 0.4, by = 0.02),
  mdrqa_radius = 0.09,
  num_samples = 100,
  sample_length = 500,
  sample_skip = 100,
  csv_file_name = "lorenz_harmonic_mdrqa_comparison.csv"
)

#
# Lorenz 96 & Harmonic oscillator
#

compute_lorenz_96_harmonic_mdrqa_comparison(
  rr = 0.5,
  coupling = seq(0, 0.4, by = 0.02),
  mdrqa_radius = 0.09,
  num_samples = 100,
  sample_length = 500,
  sample_skip = 100,
  csv_file_name = "lorenz_96_harmonic_mdrqa_comparison.csv"
)


