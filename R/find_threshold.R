## ---------------------------
##
## File name: find_threshold.R
##
## Purpose: Fast algorithm to find radius to give fixed recurrence rate
##
## Author: Sebastian Wallot & Dan Moenster
##
## Date Created: 30 May 2025
##
## ---------------------------
##
## Notes:
##
## The find_threshold function uses the distance matrix to find
## a distance (threshold) that gives a certain recurrence rate.
##
## This is used in mvjrqa() with the option setrec = TRUE.
##
## ---------------------------

find_threshold <- function(ts, target_pct, tolerance = 0.5,
                           rescale = 0, normalize = 0,
                           delay = 1, embed = 1) {
  # Ensure ts is a data frame
  ts <- as.data.frame(ts)
  # Time-delay embedding (if embed > 1)
  if (embed > 1) {
    n_rows <- nrow(ts)
    max_shift <- delay * (embed - 1)
    # Pre-check: make sure we have enough rows
    if (n_rows <= max_shift) {
      stop("Not enough rows for time-delay embedding with current delay and
           embed values.")
    }
    # Initialize embedded data with the "base" (earliest aligned) segment
    embedded_data <- ts[(max_shift + 1):n_rows, ]
    # Append each delayed version
    for (i in 1:(embed - 1)) {
      shift <- i * delay
      # Delayed part starts at (max_shift + 1 - shift) so columns align at end
      delayed_part <- ts[(max_shift + 1 - shift):(n_rows - shift), ]
      embedded_data <- cbind(embedded_data, delayed_part)
    }
    ts <- embedded_data
  }
  # Re-scale input columns if needed
  ts <- switch(
    as.character(normalize),
    "0" = ts,
    "1" = as.data.frame(lapply(ts, function(col) {
      rng <- range(col, na.rm = TRUE)
      if (rng[1] == rng[2]) {
        rep(0, length(col))
      } else {
        (col - rng[1]) / (rng[2] - rng[1])
      }
    })),
    "2" = as.data.frame(scale(ts)),
    stop("Invalid normalize option: must be 0, 1, or 2")
  )
  # Compute pairwise distances
  dists <- dist(as.matrix(ts))
  rm(ts)
  # Re-scale distances if needed
  dists <- switch(as.character(rescale),
                  `0` = dists,
                  `1` = dists / mean(dists),
                  `2` = dists / max(dists),
                  stop("Invalid rescale option: must be 0, 1, or 2")
  )
  # Define objective function for root finding
  f <- function(thresh) mean(dists < thresh) * 100 - target_pct
  # Ensure dists has a non-zero range before applying uniroot
  if (length(unique(dists)) < 2) {
    stop("Distance values are too uniform for threshold estimation.")
  }
  # Solve for threshold using root-finding
  result <- uniroot(
    f,
    c(0.05, 50),
    lower = 0.00001,
    maxiter = 25,
    extendInt = "yes",
    trace = TRUE)
  rad <- result$root
  # Check whether rad is bigger than 0
  if (rad <= 0) {
    rad <-  min(dists[dists > 0])
  }
  # Check if the recurrence rate is within the given tolerance
  # and issue a warning if it is not.
  if (abs(result$f.root) > tolerance) {
    rr_rounding <- max(0, 1 - floor(log10(tolerance)))
    rr_achieved <- round(result$f.root + target_pct, rr_rounding)
    warning("Unable to set radius to give RR within tolerance\n",
            "  Target RR = ", target_pct, "\n",
            "  Achieved RR = ", rr_achieved, "\n",
            "  Radius r = ", rad, "\n",
            "  This is a warning, please check your results.",
            call. = TRUE
            )
  }
  return(rad)
}
