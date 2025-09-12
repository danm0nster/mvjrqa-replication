## ---------------------------
##
## File name: utils.R
##
## Purpose:  Utility functions
##
## Author: Dan Moenster & Sebastian Wallot 
##
## Date Created: Apr 28 2022
##
## ---------------------------
##
## Notes:
##
## This file contains various functions used by mvjrqa()
## as well as function used to compute bootstrapped
## confidence intervals of the results.
##
## ---------------------------

library(ggplot2)
library(dplyr)
library(crqa)
library(Matrix)

# Compute recurrence rate from recurrence plot matrix.
# Needed, e.g., for the joint recurrence plot.
recurrence_rate <- function(rp) {
  return(sum(rp) / length(rp) * 100)
}

# Function to infer RQA type from time series dimensions
infer_rqa_type <- function(time_series) {
  # If time series is a vector dim(time_series) will return NULL.
  # Solution: cast to matrix first to check.
  if (dim(as.matrix(time_series))[2] == 1) {
    rqa_type <- "rqa"
  } else {
    rqa_type <- "mdcrqa"
  }
  rqa_type
}

# Function to compute RQA of one time series
rqa <- function(ts, delay, embed, radius) {
  rqa_ts <-  crqa(ts, ts,
                  delay, embed, radius,
                  rescale     = 0,           # 0 (none) | 1 (mean distance) | 2 (max distance) | 3 (min distance) | 4 (euclidean distance)
                  normalize   = 2,           # 0 (none) | 1 (unit interval) | 2 (z-score)
                  mindiagline = 2,           # min recurrent points for diagonal lines
                  minvertline = 2,           # min recurrent points for vertical lines
                  tw          = 1,           # Theiler window: 0 (incl LOI) | 1 (excl LOI)
                  whiteline   = FALSE,       # calculate empty vertical lines or not  # TRUE gives error in crqa::tt() internal helper
                  recpt       = FALSE,       # calculate measures directly from RP or not
                  side        = "both",      # base measures on "upper" or "lower" triangle, or "both"
                  method      = "rqa",       # "rqa" | "crqa" | "mdcrqa"
                  metric      = "euclidean", # distance metric: "euclidean" | "maximum" | "minkowski" | ...
                  datatype    = "continuous")
  
  return(rqa_ts)
}

# Function to compute MdRQA of one time series
mdrqa <- function(ts, delay, embed, radius) {
  rqa_ts <-  crqa(ts, ts,
                    delay, embed, radius,
                    rescale     = 0,           # 0 (none) | 1 (mean distance) | 2 (max distance) | 3 (min distance) | 4 (euclidean distance)
                    normalize   = 2,           # 0 (none) | 1 (unit interval) | 2 (z-score)
                    mindiagline = 2,           # min recurrent points for diagonal lines
                    minvertline = 2,           # min recurrent points for vertical lines
                    tw          = 1,           # Theiler window: 0 (incl LOI) | 1 (excl LOI)
                    whiteline   = FALSE,       # calculate empty vertical lines or not  # TRUE gives error in crqa::tt() internal helper
                    recpt       = FALSE,       # calculate measures directly from RP or not
                    side        = "both",      # base measures on "upper" or "lower" triangle, or "both"
                    method      = "mdcrqa",    # "rqa" | "crqa" | "mdcrqa"
                    metric      = "euclidean", # distance metric: "euclidean" | "maximum" | "minkowski" | ...
                    datatype    = "continuous")
  
  return(rqa_ts)
}

# Function to compute RQA measures for RP1, RP2 and the resulting joint RP
compute_jrqa_measures <- function(rp1, rp2,
                                  mindiagline = 2,
                                  minvertline = 2) {
   rp_list <- list(
    as.matrix(rp1),
    as.matrix(rp2),
    as.matrix(rp1) * as.matrix(rp2)
  )
  # Create results list
  results_list <- list()
  # Get recurrence measures for RPs and JRP:
  for (i in 1:3) {
    # store dimensions of the matrix
    crp <- rp_list[[i]]
    nrw <- nrow(crp)
    ncl <- ncol(crp)
    crp[is.na(crp)] <- 0
    # Computing the line counts
    numrecurs <- length(which(crp == TRUE))
    # check here if there any recurrences at all
    if (numrecurs > 0) { ## there is nothing
      RR <- sum(crp) / length(crp) * 100
      # calculate diagonal and vertical line distributions
      diagLine <- split(crp, row(crp) - col(crp))
      diaglines <- sort(
        as.numeric(unlist(lapply(
          diagLine,
          function(x) {
            runs <- rle(x)
            lx <- runs$lengths[runs$values == 1]
            return(lx)
          }
        ))),
        decreasing = TRUE
      )
      # delete line counts less than the minimum diagonal
      dcrit <- which(diaglines < mindiagline)
      if (length(dcrit) > 0) {
        diaglines <- diaglines[-dcrit]
      }
      if (length(diaglines) != 0) {
        DNRLINES <- length(diaglines) # extract the number of diagonal lines
        MDL <- max(diaglines) # extract the max length of diagonal lines
        ADL <- mean(diaglines) # compute the mean of the diagonal lines
        tabled <- as.data.frame(table(diaglines))
        total <- sum(tabled$Freq)
        p <- tabled$Freq / total
        # remove zero probability
        del <- which(p == 0)
        if (length(del) > 0) {
          p <- p[-del]
        }
        # compute entropy of diagonal line distribution
        DENTR <- -sum(p * log(p))
        # compute ratio of diagonally connected lines to recurrence points
        DET <- sum(diaglines) / numrecurs * 100
      } else {
        DNRLINES <- 0
        MDL <- 0
        ADL <- 0
        DENTR <- 0
        DET <- 0
      }
      vertLine <- split(crp, col(crp))
      vertlines <- sort(
        as.numeric(unlist(lapply(
          vertLine,
          function(x) {
            runs <- rle(x)
            lx <- runs$lengths[runs$values == 1]
            ## we could return also the black lines by just changing the above
            ## parameter to 0
            return(lx)
          }
        ))),
        decreasing = TRUE
      )
      # delete line counts less than the minimum vertical.
      vcrit <- which(vertlines < minvertline)
      if (length(vcrit) > 0) {
        vertlines <- vertlines[-vcrit]
      }
      if (length(vertlines) != 0) {
        VNRLINES <- length(vertlines) # extract the number of vertical lines
        MVL <- max(vertlines) # extract the max length of vertical lines
        AVL <- mean(vertlines) # compute the mean of the vertical lines
        tabled <- as.data.frame(table(vertlines))
        total <- sum(tabled$Freq)
        p <- tabled$Freq / total
        # remove zero probability
        del <- which(p == 0)
        if (length(del) > 0) {
          p <- p[-del]
        }
        # compute entropy of vertical line distribution
        VENTR <- -sum(p * log(p))
        # compute ratio of vertically connected lines to recurrence points
        LAM <- sum(vertlines) / numrecurs * 100
      } else {
        VNRLINES <- 0
        MVL <- 0
        AVL <- 0
        VENTR <- 0
        LAM <- 0
      }
    } else { # in case you find no recurrence
      RR <- 0 # DAN, 21 Jul 2024: This was missing. Crucial bug!
      DNRLINES <- 0
      MDL <- 0
      ADL <- 0
      DET <- NA
      DENTR <- NA
      VNRLINES <- 0
      MVL <- 0
      AVL <- 0
      LAM <- NA
      VENTR <- NA
      RP <- NA
    }
    # make crp a sparse matrix
    rp_list[[i]] <- Matrix(crp, sparse = TRUE)
    results_list[[i]] <- list(
      RR = RR,
      DET = DET,
      DNRLINES = DNRLINES,
      MDL = MDL,
      ADL = ADL,
      DENTR = DENTR,
      LAM = LAM,
      AVL = AVL,
      MVL = MVL,
      VNRLINES = VNRLINES,
      VENTR = VENTR
    )
  }
  results_list$rp_list <- rp_list
  results_list
}

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
      return(t.test(x, mu = mean(x, na.rm = na.rm), conf.level = level)$conf.int)
    } else if (method == "bootstrap") {
      bs_result <- replicate(n_samples, mean(sample(x, replace = TRUE), na.rm = na.rm))
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