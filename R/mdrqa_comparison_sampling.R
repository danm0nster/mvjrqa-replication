## ---------------------------
##
## File name: mdrqa_comparison_sampling.R
##
## Purpose:  Functions for generating multiple samples of model time series
##           used to compare MdRQA to MvJRQA.
##
## Author: Dan Moenster & Sebastian Wallot
##
## Date Created: Jun 3 2025
##
## ---------------------------
##
## Notes:
##
## Functions for generating sample data from models with random initial
## conditions.
##
## ---------------------------

# load libraries
library(crqa)
library(readr)

# source functions
source("R/mvjrqa.R")
source("R/utils.R")
source("R/lorenz_harmonic.R")
source("R/coupled_logistic_map.R")
source("R/lorenz96.R")

#
# Sampler for linear stochastic system
#
compute_linear_stochastic_mdrqa_comparison <- function(
    coupling = 1,
    rr = 20,
    mdrqa_radius = 0.6,
    num_samples = 1, # Set to 100 in script
    sample_length = 100, # Should be 500 like for the other systems!?
    csv_file_name = NULL
) {
  # Check if file name is provided
  if (!is.null(csv_file_name)) {
    # A file name has been provided
    full_file_name <- paste0("data/", csv_file_name)
  } else {
    stop("No file name provided.")
  }
  # loop through computation and save results
  k_values <- coupling
  no_samples <- num_samples
  sample_size <- sample_length
  delay_1 <- 1
  embed_1 <- 1
  radius_1 <- 0.3
  delay_2 <- 1
  embed_2 <- 1
  radius_2 <- 0.3
  delay_3 <- 1
  embed_3 <- 1
  radius_3 <- mdrqa_radius
  coupling <- NULL
  index_simu <- NULL
  rqa_RR1 <- NULL
  rqa_RR2 <- NULL
  rqa_RR12 <- NULL
  jrqa_RR <- NULL
  jrqa_DET <- NULL
  jrqa_DENTR <- NULL
  mdrqa_RR <- NULL
  mdrqa_DET <- NULL
  mdrqa_DENTR <- NULL
  k <- 0

  for (j in 1:length(k_values)) {
    for (i in 1:no_samples) {
      k <- k + 1
      # create data
      epsilon_4 <- rnorm(sample_size)
      x_1 <- rnorm(sample_size)
      y_1 <- k_values[j] * x_1 + epsilon_4 + rnorm(sample_size)
      y_2 <- k_values[j] * x_1 + epsilon_4 + rnorm(sample_size)
      Y <- as.matrix(cbind(y_1, y_2))
      ALL <- as.matrix(cbind(y_1, y_2, x_1))
      # run mvjrqa
      res_mvjrqa <- mvjrqa(as.matrix(ALL[, 1:2]),
                           as.matrix(ALL[, 3]),
                           delay_1 = 1,
                           embed_1 = 1,
                           radius_1 = .1,
                           delay_2 = 1,
                           embed_2 = 1,
                           radius_2 = .1,
                           rescale = 0,
                           normalize = 2,
                           mindiagline = 2,
                           minvertline = 2,
                           tw = 1,
                           whiteline = TRUE,
                           recpt = FALSE,
                           side = "both",
                           metric = "euclidean",
                           datatype = "continuous",
                           setrec = TRUE,
                           targetrec = rr)
      # run mdrqa
      res_mdrqa <- crqa(ALL, ALL,
                       delay_3, embed_3, radius_3,
                       rescale     = 0,
                       normalize   = 2,
                       mindiagline = 2,
                       minvertline = 2,
                       tw          = 1,
                       whiteline   = FALSE,
                       recpt       = FALSE,
                       side        = "both",
                       method      = "mdcrqa",
                       metric      = "euclidean",
                       datatype    = "continuous")
      ## collect data
      # get beta value
      coupling[k] <- k_values[j]
      # get simulation trial
      index_simu[k] <- k
      # get RR of first set of time series
      rqa_RR1[k] <- res_mvjrqa[[1]]$RR
      # get RR of second set of time series
      rqa_RR2[k] <- res_mvjrqa[[2]]$RR
      # get average individual RR
      rqa_RR12[k] <- (res_mvjrqa[[1]]$RR + res_mvjrqa[[2]]$RR) / 2
      # get joint RR
      jrqa_RR[k] <- res_mvjrqa[[3]]$RR
      # get joint DET
      jrqa_DET[k] <- res_mvjrqa[[3]]$DET
      # get joint ENTR
      jrqa_DENTR[k] <- res_mvjrqa[[3]]$DENTR
      # get mdrqa RR
      mdrqa_RR[k] <- res_mdrqa$RR
      # get mdrqa DET
      mdrqa_DET[k] <- res_mdrqa$DET
      # get mdrqa ENTR
      mdrqa_DENTR[k] <- res_mdrqa$ENTR
    }
  }
  # collect data
  joint_rr <- data.frame(coupling, index_simu, rqa_RR1, rqa_RR2, rqa_RR12,
                         jrqa_RR, jrqa_DET, jrqa_DENTR, mdrqa_RR, mdrqa_DET,
                         mdrqa_DENTR)
  # Write results to file
  create_dir_if_not_present("data")
  write_csv(joint_rr, file = full_file_name)
}


#
# Coupled logistic maps with external driver
#

compute_logistic_system_mdrqa_comparison <- function(coupling = 1,
                                                     rr = 10,
                                                     mdrqa_radius = 0.5,
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
                                  whiteline = TRUE,
                                  recpt = FALSE,
                                  side = "both",
                                  metric = "euclidean",
                                  datatype = "continuous"
        )
        rqa_system_1 <- as.data.frame(logistic_mvjrqa[[1]]) |>
          mutate(sample_id = sample_number,
                 method = "MvJRQA",
                 system = "1",
                 coupling = cc,
                 RR_target = rr_target)
        rqa_system_2 <- as.data.frame(logistic_mvjrqa[[2]]) |>
          mutate(sample_id = sample_number,
                 method = "MvJRQA",
                 system = "2",
                 coupling = cc,
                 RR_target = rr_target)
        rqa_joint <- as.data.frame(logistic_mvjrqa[[3]]) |>
          mutate(sample_id = sample_number,
                 method = "MvJRQA",
                 system = "j",
                 coupling = cc,
                 RR_target = rr_target)
        #
        # Calculate MdRQA
        #
        full_ts <- as.matrix(model_data)
        res_mdrqa <- crqa(full_ts, full_ts,
                         delay, embed, mdrqa_radius,
                         rescale     = 0,
                         normalize   = 2,
                         mindiagline = 2,
                         minvertline = 2,
                         tw          = 1,
                         whiteline   = TRUE,
                         recpt       = FALSE,
                         side        = "both",
                         method      = "mdcrqa",
                         metric      = "euclidean",
                         datatype    = "continuous")
        # Compute them using our function with the same rp twice
        mdrqa_results <- as.data.frame(
          compute_jrqa_measures(res_mdrqa$RP,
                                res_mdrqa$RP)[[3]]
        ) |>
          mutate(sample_id = sample_number,
                 method = "MdRQA",
                 system = "f",
                 coupling = cc,
                 radius = mdrqa_radius,
                 RR_target = NA)
        # Add rows to joint_rr data frame
        joint_rr <- bind_rows(joint_rr,
                              rqa_system_1,
                              rqa_system_2,
                              rqa_joint,
                              mdrqa_results)
      }
    }
  }
  create_dir_if_not_present("data")
  write_csv(joint_rr, file = full_file_name)
}



#
# Lorenz-Harmonic oscillator
#

compute_lorenz_harmonic_mdrqa_comparison <- function(coupling = 1,
                                                     rr = 1,
                                                     mdrqa_radius = 0.1,
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
  # Check if file name is provided
  if(!is.null(csv_file_name)) {
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
        #
        # Calculate mvjrqa with fixed sub-system RR
        #
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
                 method = "MvJRQA",
                 system = "1",
                 coupling = cc,
                 RR_target = rr_target)
        rqa_system_2 <- as.data.frame(nonlinear_mvjrqa[[2]]) |>
          mutate(sample_id = sample_number,
                 method = "MvJRQA",
                 system = "2",
                 coupling = cc,
                 RR_target = rr_target)
        rqa_joint <- as.data.frame(nonlinear_mvjrqa[[3]]) |>
          mutate(sample_id = sample_number,
                 method = "MvJRQA",
                 system = "j",
                 coupling = cc,
                 RR_target = rr_target)
      
        #
        # Calculate MdRQA
        #
        full_ts <- as.matrix(model_data)
        res_mdrqa <- crqa(full_ts, full_ts,
                         delay, embed, mdrqa_radius,
                         rescale     = 0,
                         normalize   = 2,
                         mindiagline = 2,
                         minvertline = 2,
                         tw          = 1,
                         whiteline   = TRUE,
                         recpt       = FALSE,
                         side        = "both",
                         method      = "mdcrqa",
                         metric      = "euclidean",
                         datatype    = "continuous")
        # Compute them using our function with the same rp twice
        mdrqa_results <- as.data.frame(
          compute_jrqa_measures(res_mdrqa$RP,
                                res_mdrqa$RP)[[3]]
        ) |>
          mutate(sample_id = sample_number,
                 method = "MdRQA",
                 system = "f",
                 coupling = cc,
                 radius = mdrqa_radius,
                 RR_target = NA)
        # Get the RQA measured directly
        mdrqa_results_2 <- data.frame(
          RR = res_mdrqa$RR,
          DET = res_mdrqa$DET,
          DNRLINES = res_mdrqa$NRLINE,
          MDL = res_mdrqa$maxL,
          ADL = res_mdrqa$L,
          DENTR = res_mdrqa$ENTR,
          LAM = res_mdrqa$LAM,
          AVL = NA, # Not computed by crqa()?
          MVL = NA,
          VNRLINES = NA,
          VENTR = NA
        ) |>
          mutate(sample_id = sample_number,
                 method = "MdRQA_crqa",
                 system = "f",
                 coupling = cc,
                 radius = mdrqa_radius,
                 RR_target = NA)
        # Add rows to joint_rr data frame
        joint_rr <- bind_rows(joint_rr,
                              rqa_system_1,
                              rqa_system_2,
                              rqa_joint,
                              mdrqa_results,
                              mdrqa_results_2)
      }
    }
  }
  create_dir_if_not_present("data")
  write_csv(joint_rr, file = full_file_name)
}

#
# Lorenz 96 & Harmonic oscillator
#

compute_lorenz_96_harmonic_mdrqa_comparison <- function(coupling = 1,
                                                        rr = 1,
                                                        n_var = 5,
                                                        mdrqa_radius = 0.1,
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
  # Check if file name is provided
  if (!is.null(csv_file_name)) {
    # A file name has been provided
    full_file_name <- paste0("data/", csv_file_name)
  } else {
    stop("No file name provided.")
  }
  joint_rr <- data.frame()
  sample_number <- 0
  # Set up random initial states for each sample
  x_init <- matrix(nrow = num_samples, ncol = n_var)
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
                                    coupling = cc,
                                    initial_state = list(x = x_init[s, ],
                                                         u = u_init[s],
                                                         v = v_init[s]))
      # Extract the two time series from model_data
      ts_1 <- extract_lorenz96(model_data)
      ts_2 <- extract_oscillator(model_data)
      for (rr_target in rr_values) {
        sample_number <- sample_number + 1
        #
        # Calculate mvjrqa with fixed sub-system RR
        #
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
                 method = "MvJRQA",
                 system = "1",
                 coupling = cc,
                 RR_target = rr_target)
        rqa_system_2 <- as.data.frame(nonlinear_mvjrqa[[2]]) |>
          mutate(sample_id = sample_number,
                 method = "MvJRQA",
                 system = "2",
                 coupling = cc,
                 RR_target = rr_target)
        rqa_joint <- as.data.frame(nonlinear_mvjrqa[[3]]) |>
          mutate(sample_id = sample_number,
                 method = "MvJRQA",
                 system = "j",
                 coupling = cc,
                 RR_target = rr_target)
      
        #
        # Calculate MdRQA
        #
        full_ts <- as.matrix(model_data)
        res_mdrqa <- crqa(full_ts, full_ts,
                         delay, embed, mdrqa_radius,
                         rescale     = 0,
                         normalize   = 2,
                         mindiagline = 2,
                         minvertline = 2,
                         tw          = 1,
                         whiteline   = TRUE,
                         recpt       = FALSE,
                         side        = "both",
                         method      = "mdcrqa",
                         metric      = "euclidean",
                         datatype    = "continuous")
        # Compute them using our function with the same rp twice
        mdrqa_results <- as.data.frame(
          compute_jrqa_measures(res_mdrqa$RP,
                                res_mdrqa$RP)[[3]]
        ) |>
          mutate(sample_id = sample_number,
                 method = "MdRQA",
                 system = "f",
                 coupling = cc,
                 radius = mdrqa_radius,
                 RR_target = NA)
        # Get the RQA measured directly
        mdrqa_results_2 <- data.frame(
          RR = res_mdrqa$RR,
          DET = res_mdrqa$DET,
          DNRLINES = res_mdrqa$NRLINE,
          MDL = res_mdrqa$maxL,
          ADL = res_mdrqa$L,
          DENTR = res_mdrqa$ENTR,
          LAM = res_mdrqa$LAM,
          AVL = NA, # Not computed by crqa()?
          MVL = NA,
          VNRLINES = NA,
          VENTR = NA
        ) |>
          mutate(sample_id = sample_number,
                 method = "MdRQA_crqa",
                 system = "f",
                 coupling = cc,
                 radius = mdrqa_radius,
                 RR_target = NA)
        # Add rows to joint_rr data frame
        joint_rr <- bind_rows(joint_rr,
                              rqa_system_1,
                              rqa_system_2,
                              rqa_joint,
                              mdrqa_results,
                              mdrqa_results_2)
      }
    }
  }
  create_dir_if_not_present("data")
  write_csv(joint_rr, file = full_file_name)
}
