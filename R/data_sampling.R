## ---------------------------
##
## File name: data_sampling.R
##
## Purpose:  Functions for generating multiple samples of model time series.
##
## Author: Dan Moenster & Sebastian Wallot
##
## Date Created: May 14 2025
##
## ---------------------------
##
## Notes:
##
## Functions for generating sample data from models with random initial
## conditions.
##
## ---------------------------

library(readr)
library(dplyr)

source("R/lorenz_harmonic.R")
source("R/lorenz96.R")
source("R/coupled_logistic_map.R")
source("R/mvjrqa.R")
source("R/utils.R")

compute_lorenz_harmonic_data <- function(coupling = 1,
                                         rr = 10,
                                         num_samples = 1,
                                         sample_length = 500,
                                         sample_skip = 100,
                                         csv_file_name = NULL) {
  rr_values <- rr
  coupling_values <- coupling
  # None of the time series need to be further embedded,
  # so delay and embed are set to 1.
  delay <- 1
  embed <- 1
  # Check if file name is provided and if data should be computed
  if (!is.null(csv_file_name)) {
    # A file name has been provided
    full_file_name <- paste0("data/", csv_file_name)
  } else {
    stop("No file name provided.")
  }

  joint_rr <- data.frame()
  sample_number <- 0
  # Set up random initial states for each sample
  x_init <- runif(num_samples, min = -10, max = 10)
  y_init <- runif(num_samples, min = -10, max = 10)
  z_init <- runif(num_samples, min = 0, max = 40)
  u_init <- runif(num_samples, min = -10, max = 10)
  v_init <- runif(num_samples, min = -20, max = 20)
  for (s in 1:num_samples) {
    for (cc in coupling_values) {
      # Run the model to produce data
      model_data <- lorenz_harmonic(n = sample_length,
                                    skip = sample_skip,
                                    coupling = cc,
                                    initial_state = c(x = x_init[s],
                                                      y = y_init[s],
                                                      z = z_init[s],
                                                      u = u_init[s],
                                                      v = v_init[s]))
      # Extract the two time series from model_data
      ts_1 <- extract_lorenz(model_data)
      ts_2 <- extract_oscillator(model_data)
      for (rr_target in rr_values) {
        sample_number <- sample_number + 1
        nonlinear_mvjrqa <- mvjrqa(ts_1,
                                   ts_2,
                                   delay_1 = delay,
                                   embed_1 = embed,
                                   radius_1 = NA,
                                   delay_2 = delay,
                                   embed_2 = embed,
                                   radius_2 = NA,
                                   setrec = TRUE,
                                   targetrec = rr_target)
        rqa_system_1 <- as.data.frame(nonlinear_mvjrqa[[1]]) |>
          mutate(sample_id = sample_number,
                 system = "1",
                 coupling = cc,
                 RR_target = rr_target)
        rqa_system_2 <- as.data.frame(nonlinear_mvjrqa[[2]]) |>
          mutate(sample_id = sample_number,
                 system = "2",
                 coupling = cc,
                 RR_target = rr_target)
        rqa_joint <- as.data.frame(nonlinear_mvjrqa[[3]]) |>
          mutate(sample_id = sample_number,
                 system = "j",
                 coupling = cc,
                 RR_target = rr_target)
        # Add rows to joint_rr data frame
        joint_rr <- bind_rows(joint_rr,
                              rqa_system_1,
                              rqa_system_2,
                              rqa_joint)
      }
    }
  }
  create_dir_if_not_present("data")
  write_csv(joint_rr, file = full_file_name)
}

compute_lorenz96_harmonic_data <- function(coupling = 1,
                                           rr = 10,
                                           n_var = 5,
                                           num_samples = 1,
                                           sample_length = 500,
                                           sample_skip = 100,
                                           csv_file_name = NULL) {
  rr_values <- rr
  coupling_values <- coupling
  # None of the time series need to be further embedded,
  # so delay and embed are set to 1.
  delay <- 1
  embed <- 1
  # Check if file name is provided and if data should be computed
  if (!is.null(csv_file_name)) {
    # A file name has been provided
    full_file_name <- paste0("data/", csv_file_name)
  } else {
    stop("No file name provided.")
  }

  joint_rr <- data.frame()
  sample_number <- 0
  x_init <- matrix(nrow = num_samples, ncol = n_var)
  # Set up random initial states for each sample
  for (x_col in 1:n_var) {
    x_init[, x_col] <- rnorm(num_samples, mean = 1, sd = 0.01)
  }
  u_init <- runif(num_samples, min = -10, max = 10)
  v_init <- runif(num_samples, min = -20, max = 20)
  for (s in 1:num_samples) {
    for (cc in coupling_values) {
      # Run the model to produce data
      model_data <- lorenz96_harmonic(n = sample_length,
                                      skip = sample_skip,
                                      n_var = n_var,
                                      coupling = cc,
                                      initial_state = list(x = x_init[s, ],
                                                           u = u_init[s],
                                                           v = v_init[s]))
      # Extract the two time series from model_data
      ts_1 <- extract_lorenz96(model_data)
      ts_2 <- extract_oscillator(model_data)
      for (rr_target in rr_values) {
        sample_number <- sample_number + 1
        nonlinear_mvjrqa <- mvjrqa(ts_1,
                                   ts_2,
                                   delay_1 = delay,
                                   embed_1 = embed,
                                   radius_1 = NA,
                                   delay_2 = delay,
                                   embed_2 = embed,
                                   radius_2 = NA,
                                   setrec = TRUE,
                                   targetrec = rr_target)
        rqa_system_1 <- as.data.frame(nonlinear_mvjrqa[[1]]) |>
          mutate(sample_id = sample_number,
                 system = "1",
                 coupling = cc,
                 RR_target = rr_target)
        rqa_system_2 <- as.data.frame(nonlinear_mvjrqa[[2]]) |>
          mutate(sample_id = sample_number,
                 system = "2",
                 coupling = cc,
                 RR_target = rr_target)
        rqa_joint <- as.data.frame(nonlinear_mvjrqa[[3]]) |>
          mutate(sample_id = sample_number,
                 system = "j",
                 coupling = cc,
                 RR_target = rr_target)
        # Add rows to joint_rr data frame
        joint_rr <- bind_rows(joint_rr,
                              rqa_system_1,
                              rqa_system_2,
                              rqa_joint)
      }
    }
  }
  create_dir_if_not_present("data")
  write_csv(joint_rr, file = full_file_name)
}


compute_linear_stochastic_system_data <- function(coupling = 1,
                                                  rr = 10,
                                                  num_samples = 100,
                                                  sample_length = 100,
                                                  csv_file_name = NULL) {
  b_1 <- coupling
  rr_values <- rr
  no_samples <- num_samples
  sample_size <- sample_length
  delay_1 <- 1
  embed_1 <- 1
  delay_2 <- 1
  embed_2 <- 1

  if (!is.null(csv_file_name)) {
    # A file name has been provided
    full_file_name <- paste0("data/", csv_file_name)
  } else {
    # No file name provided
    stop("No file name provided.")
  }
  joint_rr <- data.frame()
  sample_number <- 0
  for (j in 1:length(b_1)) {
    for (i in 1:no_samples) {
      # create data
      epsilon_4 <- rnorm(sample_size)
      x_1 <- rnorm(sample_size)
      y_1 <- b_1[j] * x_1 + epsilon_4 + rnorm(sample_size)
      y_2 <- b_1[j] * x_1 + epsilon_4 + rnorm(sample_size)
      Y <- as.matrix(cbind(y_1, y_2))
      for (rr_target in rr_values) {
        sample_number <- sample_number + 1
        linear_mvjrqa <- mvjrqa(x_1, Y,
                                delay_1 = delay_1,
                                embed_1 = embed_1,
                                radius_1 = NA,
                                delay_2 = delay_2,
                                embed_2 = embed_2,
                                radius_2 = NA,
                                datatype = "continuous",
                                setrec = TRUE,
                                targetrec = rr_target)
        # Add columns for system 1, 2 and joint RP
        rqa_system_1 <- as.data.frame(linear_mvjrqa[[1]]) |>
          mutate(sample_id = sample_number,
                 system = "1",
                 coupling = b_1[j],
                 RR_target = rr_target)
        rqa_system_2 <- as.data.frame(linear_mvjrqa[[2]]) |>
          mutate(sample_id = sample_number,
                 system = "2",
                 coupling = b_1[j],
                 RR_target = rr_target)
        rqa_joint <- as.data.frame(linear_mvjrqa[[3]]) |>
          mutate(sample_id = sample_number,
                 system = "j",
                 coupling = b_1[j],
                 RR_target = rr_target)
      
        joint_rr <- bind_rows(joint_rr,
                              rqa_system_1,
                              rqa_system_2,
                              rqa_joint)
      }
    }
  }   
  create_dir_if_not_present("data")
  write_csv(joint_rr, file = full_file_name)
}


###
###  Function to compute data for the driven coupled logistic map
###
compute_logistic_system_data <- function(coupling = 1,
                                         rr = 10,
                                         rx = 3.58,
                                         ry = 3.6,
                                         byx = 0.1,
                                         bxy = 0,
                                         num_samples = 1,
                                         sample_length = 500,
                                         sample_skip = 300,
                                         csv_file_name = NULL,
                                         DEBUG = FALSE) {
  rr_values <- rr
  coupling_values <- coupling
  # Set delay and embed as in Mønster et al. (2017):
  #   Causal inference from noisy time-series data —
  #   Testing the Convergent Cross-Mapping algorithm
  #   in the presence of noise and external influence
  delay <- 1
  embed <- 2
  # Both the driver and the logistic maps will be embedded
  # Note: This means we have a 2-dimensional and a 4-dimensional system in
  # phase space.
  # Check if file name is provided and if data should be computed
  if (!is.null(csv_file_name)) {
    # A file name has been provided
    full_file_name <- paste0("data/", csv_file_name)
  } else {
    # No file name provided
    stop("No file name provided.")
  }

  joint_rr <- data.frame()
  sample_number <- 0
  # Set up random initial states for each sample
  x_init <- runif(num_samples, min = 0.2, max = 0.8)
  y_init <- runif(num_samples, min = 0.2, max = 0.8)
  phi_init <- runif(num_samples, min = 0, max = 2 * pi)
  for (s in 1:num_samples) {
    for (rr_target in rr_values) {
      for (cc in coupling_values) {
        # Run the model to produce data
        if (DEBUG) {
          print("Arguments:")
          print("")
          print(paste("x0 = ", x_init[s]))
          print(paste("y0 = ", y_init[s]))
          print(paste("rx =",  rx))
          print(paste("ry = ", ry))
          print(paste("bxy = ", bxy))
          print(paste("byx = ", byx))
          print(paste("etax = ", cc))
          print(paste("etay = ", cc))
          print(paste("p = ", 30))
          print(paste("phi = ", phi_init[s]))
          print(paste("N =", sample_length))
          print(paste("N_skip = ", sample_skip))
        }
        model_data <- coupled_logistic_map(
          x0 = x_init[s],
          y0 = y_init[s],
          rx = rx,
          ry = ry,
          bxy = bxy,
          byx = byx,
          etax = cc,
          etay = cc,
          p = 30,
          phi = phi_init[s],
          N = sample_length,
          N_skip = sample_skip
        )
        # Extract the two time series from model_data
        ts_1 <- as.matrix(extract_H(model_data))
        ts_2 <- matrix(
          c(
            extract_X(model_data),
            extract_Y(model_data)
          ),
          ncol = 2, byrow = FALSE
        )
        sample_number <- sample_number + 1
        logistic_mvjrqa <- mvjrqa(ts_1,
                                  ts_2,
                                  delay_1 = delay,
                                  embed_1 = embed,
                                  radius_1 = NA,
                                  delay_2 = delay,
                                  embed_2 = embed,
                                  radius_2 = NA,
                                  setrec = TRUE,
                                  targetrec = rr_target,
                                  rescale = 0,
                                  normalize = 2,
                                  mindiagline = 2,
                                  minvertline = 2,
                                  tw = 1,
                                  whiteline = FALSE,
                                  recpt = FALSE,
                                  side = "both",
                                  metric = "euclidean",
                                  datatype = "continuous"
        )
        rqa_system_1 <- as.data.frame(logistic_mvjrqa[[1]]) |>
          mutate(sample_id = sample_number,
                 system = "1",
                 coupling = cc,
                 RR_target = rr_target)
        rqa_system_2 <- as.data.frame(logistic_mvjrqa[[2]]) |>
          mutate(sample_id = sample_number,
                 system = "2",
                 coupling = cc,
                 RR_target = rr_target)
        rqa_joint <- as.data.frame(logistic_mvjrqa[[3]]) |>
          mutate(sample_id = sample_number,
                 system = "j",
                 coupling = cc,
                 RR_target = rr_target)
        # Add rows to joint_rr data frame
        joint_rr <- bind_rows(joint_rr,
                              rqa_system_1,
                              rqa_system_2,
                              rqa_joint)
      }
    }
  }
  create_dir_if_not_present("data")
  write_csv(joint_rr, file = full_file_name)
}
