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
## This file contains various functions used by mvjrqa().
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
                  rescale     = 0,
                  normalize   = 2,
                  mindiagline = 2,
                  minvertline = 2,
                  tw          = 1,
                  whiteline   = FALSE,
                  recpt       = FALSE,
                  side        = "both",
                  method      = "rqa",
                  metric      = "euclidean",
                  datatype    = "continuous")
  return(rqa_ts)
}

# Function to compute MdRQA of one time series
mdrqa <- function(ts, delay, embed, radius) {
  rqa_ts <-  crqa(ts, ts,
                  delay, embed, radius,
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
