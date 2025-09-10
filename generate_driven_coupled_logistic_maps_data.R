#
# Calculate data for the driven coupled logistic map system
#

source("R/data_sampling.R")


# Coupling sweep
couplings <- seq(0, 0.08, by = 0.005)
# Because system 1 is so regular (cosine), it is not possible to
# get very low RR. The lowest is around RR = 4%
rr_values <- c(4, 6, 8, 10, 16)

# Set random seed to make the script reproducible
set.seed(84435)

compute_logistic_system_data(
  coupling = couplings,
  rr = rr_values,
  num_samples = 100,
  sample_length = 500,
  sample_skip = 300,
  csv_file_name = "driven_coupled_logistic_maps_coupling_sweep.csv"
)

# Extreme couplings
extreme_couplings <- seq(0, 0.3, by = 0.02)
rr_values <- c(6, 8, 10)

# Set random seed to make the script reproducible
set.seed(69742)

compute_logistic_system_data(
  coupling = extreme_couplings,
  rr = rr_values,
  num_samples = 100,
  sample_length = 500,
  sample_skip = 300,
  csv_file_name = "driven_coupled_logistic_maps_extreme_coupling.csv"
)

# RR sweep
couplings <- c(0, 0.02, 0.04, 0.08)
# Because system 1 is so regular (cosine), it is not possible to
# get very low RR. The lowest seems to be around RR = 4%
rr_values <- c(4, 6, 10, 20, 40, 80)

# Set random seed to make the script reproducible
set.seed(59882)

compute_logistic_system_data(
  coupling = couplings,
  rr = rr_values,
  num_samples = 100,
  sample_length = 500,
  sample_skip = 300,
  csv_file_name = "driven_coupled_logistic_maps_rr_sweep.csv"
)

