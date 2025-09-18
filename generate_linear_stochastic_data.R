## ---------------------------
##
## Script name: generate_linear_stochastic_data.R
##
## Purpose of script: Calculate data for the linear stochastic system
##
## Author: Dan Moenster & Sebastian Wallot
##
## Date Created: May 25 2025
##
## ---------------------------
##
## Notes:
##
## This script will produce the files:
## data/linear_stochastic_coupling_sweep.csv
## data/linear_stochastic_extreme_coupling.csv
## data/linear_stochastic_rr_sweep.csv
##
## ---------------------------

source("R/data_sampling.R")

#
# Sweep over moderate values of the coupling.
#

rr_values <- c(1, 2, 5, 10)
coupling_values <- seq(0, 10, by = 0.5)

# Set random seed to make the script reproducible
set.seed(80885)

compute_linear_stochastic_system_data(
  coupling = coupling_values,
  rr = rr_values,
  csv_file_name = "linear_stochastic_coupling_sweep.csv"
)


#
# Extreme coupling
#
coupling_values <- seq(0, 100, by = 5)
rr_values <- c(1, 2, 5, 10)

# Set random seed to make the script reproducible
set.seed(27670)

compute_linear_stochastic_system_data(
  coupling = coupling_values,
  rr = rr_values,
  csv_file_name = "linear_stochastic_extreme_coupling.csv"
)

#
# Sweep over RR-values for several values of coupling
#

coupling_values <- c(0, 3, 6, 12) 
rr_values <- c(0.3, 0.6, 1.25, 2.5, 5, 10, 20, 40, 80)

# Set random seed to make the script reproducible
set.seed(58373)

compute_linear_stochastic_system_data(
  coupling = coupling_values,
  rr = rr_values,
  csv_file_name = "linear_stochastic_rr_sweep.csv"
)