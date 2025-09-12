## ---------------------------
##
## File name: bootstrap.R
##
## Purpose:  Utility functions
##
## Author: Dan Moenster & Sebastian Wallot
##
## Date Created: Sep 12 2025
##
## ---------------------------
##
## Notes:
##
## This file contains  functions used to compute bootstrapped
## confidence intervals of the results.
##
## The functions were moved here from utils.R
##
## ---------------------------

library(dplyr)

#
# Function to compute bootstrapped CI from a vector of sampled values
#
# The function returns a vector (ci_lo, ci_hi)
bootstrap_ci <- function(x,
                         na.rm = TRUE,
                         level = 0.95,
                         method = "bootstrap",
                         center = TRUE,
                         samples = NULL) {
  if (is.null(samples)) {
    n_samples <- length(x)
  } else if (is.numeric(samples) && samples > 1) {
    if (ceiling(samples) == floor(samples)) {
      n_samples <- samples
    } else {
      stop("samples must be an integer")
    }
  } else {
    stop("samples must be an integer > 1")
  }
  # In case all values are NA. Return NA. This occurs for the radius
  # of the joint recurrence plot which is not defined.
  if (all(is.na(x))) {
    return(NA)
  }
  if (sd(x, na.rm = na.rm) == 0) {
    # If standard deviation is zero it makes no sense to compute a CI
    return(NA)
  } else {
    if (method == "t.test") {
      return(t.test(x,
                    mu = mean(x, na.rm = na.rm),
                    conf.level = level)$conf.int
      )
    } else if (method == "bootstrap") {
      bs_result <- replicate(n_samples,
                             mean(sample(x, replace = TRUE),
                                  na.rm = na.rm))
      bs_lo <- quantile(bs_result, probs = (1 - level) / 2)
      bs_hi <- quantile(bs_result, probs = 1 - (1 - level) / 2)
      if (center) {
        # Optionally center the CI on the original sample mean
        bs_mean <- mean(bs_result)
        raw_mean <- mean(x, na.rm = na.rm)
        ci_lo <- bs_lo + raw_mean - bs_mean
        ci_hi <- bs_hi + raw_mean - bs_mean
      } else {
        ci_lo <- bs_lo
        ci_hi <- bs_hi
      }
      return(c(ci_lo, ci_hi))
    }
  }
}
#
# Function to compute the lower endpoint of the CI.
# This function is a wrapper around bootstrap_ci()
#
lower_ci <- function(x,
                     na.rm = TRUE,
                     level = 0.95,
                     method = "bootstrap",
                     center = TRUE,
                     samples = NULL) {
  ci <- bootstrap_ci(x,
                     na.rm = na.rm,
                     level = level,
                     method = method,
                     center = center,
                     samples = samples)
  return(ci[1])
}
#
# Function to compute the upper endpoint of the CI.
# This function is a wrapper around bootstrap_ci()
#
upper_ci <- function(x,
                     na.rm = TRUE,
                     level = 0.95,
                     method = "bootstrap",
                     center = TRUE,
                     samples = NULL) {
  ci <- bootstrap_ci(x,
                     na.rm = na.rm,
                     level = level,
                     method = method,
                     center = center,
                     samples = samples)
  return(ci[2])
}
#
# Function to compute standard error of the mean from a vector of values
#
standard_error <- function(x, na.rm = TRUE) {
  if (na.rm) {
    # If NA values are removed in sd() they should not be counted
    N <- sum(!is.na(x))
  } else {
    N <- length(x)
  }
  return(sd(x, na.rm = na.rm) / sqrt(N))
}

#
# Function to compute summary of recurrence quantities from a sample
# If there is a sample_id column, that is not summarised
#
compute_sample_summary <- function(rp_sample) {
  rp_summary <- rp_sample |>
    group_by(RR_target, coupling, system) |>
    summarise(across(
      .cols = !sample_id & where(is.numeric),
      .fns = list(M = mean, SD = sd, SE = standard_error,
                  CI_lo = lower_ci, CI_hi = upper_ci), na.rm = TRUE,
      .names = "{col}_{fn}"
    ), .groups = "drop")
  return(rp_summary)
}
