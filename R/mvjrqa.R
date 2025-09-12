## ---------------------------
##
## File name: mvjrqa.R
##
## Purpose: Main functions for MvJRQA analysis
##
## Author: Dan Moenster & Sebastian Wallot 
##
## Date Created: Apr 12 2024
##
## ---------------------------
##
## Notes:
##
## This file contains the main mvjrqa() function that will call
## either fixed_radius_mvjrqa() or fixed_recurrence_mvjrqa()
## depending on whether the value of setrec is FALSE or TRUE.
##
## ---------------------------
library(crqa)
library(Matrix)
# get helper functions
source("R/utils.R")
source("R/find_threshold.R")

# fixed_radius_mvjrqa() function  ===========================
#
# Function to perform mvjrqa for fixed radii. This is called by mvjrqa().
#

fixed_radius_mvjrqa <- function(ts_1,
                                ts_2,
                                delay_1,
                                embed_1,
                                radius_1,
                                delay_2,
                                embed_2,
                                radius_2,
                                rescale = 0, # 0 is "none". Allow character
                                normalize = 2, # 2 is "z-score". Allow character
                                mindiagline = 2,
                                minvertline = 2,
                                tw = 1,
                                whiteline = FALSE,  # TRUE gives error in crqa::tt() internal helper
                                recpt = FALSE,
                                side = "both",
                                metric = "euclidean",
                                datatype = "continuous") {
  # run rqa for ts_1
  #
  # Use the dimensions of time_series to infer the recurrence method to call
  # crqa() with.
  rqa_type <- infer_rqa_type(ts_1)
  if (rqa_type == "rqa") {
    # crqa() does not want matrices with only one column, so cast to vector.
    ts_1 <- as.vector(ts_1)
  }
  rqa_ts_1 <- crqa(
    ts1 = ts_1, ts2 = ts_1,
    delay_1, embed_1, radius_1,
    rescale = rescale,
    normalize = normalize,
    mindiagline = mindiagline,
    minvertline = minvertline,
    tw = tw,
    whiteline = whiteline,
    recpt = recpt,
    side = side,
    method = rqa_type,
    metric = metric,
    datatype = datatype
  )
  # run rqa for ts_2
  rqa_type <- infer_rqa_type(ts_2)
  if (rqa_type == "rqa") {
    # crqa() does not want matrices with only one column, so cast to vector.
    ts_2 <- as.vector(ts_2)
  }
  rqa_ts_2 <- crqa(
    ts1 = ts_2, ts2 = ts_2,
    delay_2, embed_2, radius_2,
    rescale = rescale,
    normalize = normalize,
    mindiagline = mindiagline,
    minvertline = minvertline,
    tw = tw,
    whiteline = whiteline,
    recpt = recpt,
    side = side,
    method = rqa_type,
    metric = metric,
    datatype = datatype
  )
  # Compute the RQA measures for RP1, RP2 and the joint RP
  results_list <- compute_jrqa_measures(rqa_ts_1$RP, rqa_ts_2$RP)
  # Add radius to results_list
  results_list[[1]]$radius <- radius_1
  results_list[[2]]$radius <- radius_2
  results_list
}

# fixed_recurrence_mvjrqa() function  ===========================
#
# Function to perform mvjrqa for fixed recurrence rate.
# This is called by mvjrqa().
#
fixed_recurrence_mvjrqa <- function(ts1,
                                    ts2,
                                    delay_1,
                                    embed_1,
                                    radius_1,
                                    delay_2,
                                    embed_2,
                                    radius_2,
                                    rescale = 0, # 0 is "none". Allow character
                                    normalize = 2, # 2 is "z-score".
                                    mindiagline = 2,
                                    minvertline = 2,
                                    tw = 1,
                                    whiteline = FALSE,  # TRUE gives error in crqa::tt() internal helper
                                    recpt = FALSE,
                                    side = "both",
                                    metric = "euclidean",
                                    datatype = "continuous",
                                    targetrec = NULL) {
  # Use find_threshold to find the radii to give fixed RR = targetrec
  radius_1_solution <- find_threshold(ts1,
                                      targetrec,
                                      rescale = rescale,
                                      normalize = normalize,
                                      delay = delay_1,
                                      embed = embed_1)
  radius_2_solution <- find_threshold(ts2,
                                      targetrec,
                                      rescale = rescale,
                                      normalize = normalize,
                                      delay = delay_2,
                                      embed = embed_2)
  # Use the dimensions of time_series to infer the recurrence method to call
  # crqa() with.
  rqa_type <- infer_rqa_type(ts1)
  if (rqa_type == "rqa") {
    # crqa() does not want matrices with only one column, so cast to vector.
    ts1 <- as.vector(ts1)
  }
  rqa_ts_1 <- crqa(ts1, ts1,
                   delay_1, embed_1, radius_1_solution,
                   rescale = rescale,
                   normalize = normalize,
                   mindiagline = mindiagline,
                   minvertline = minvertline,
                   tw = tw,
                   whiteline = whiteline,
                   recpt = recpt,
                   side = side,
                   method = rqa_type,
                   metric = metric,
                   datatype = datatype
  )
  rqa_type <- infer_rqa_type(ts2)
  if (rqa_type == "rqa") {
    # crqa() does not want matrices with only one column, so cast to vector.
    ts2 <- as.vector(ts2)
  }
  rqa_ts_2 <- crqa(ts2, ts2,
                   delay_2, embed_2, radius_2_solution,
                   rescale = rescale,
                   normalize = normalize,
                   mindiagline = mindiagline,
                   minvertline = minvertline,
                   tw = tw,
                   whiteline = whiteline,
                   recpt = recpt,
                   side = side,
                   method = rqa_type,
                   metric = metric,
                   datatype = datatype
  )
  # Compute RQA measures for RP1, RP2 and the joint RP
  if (!identical(dim(rqa_ts_1$RP), dim(rqa_ts_2$RP))) {
    plot(ts_2)
    stop("Recurrence plot of ts_1 and ts_2 have different dimensions:",
         "\n  dim(RP_1) = ", paste(dim(rqa_ts_1$RP), collapse = " x "),
         "\n  dim(RP_2) = ", paste(dim(rqa_ts_2$RP), collapse = " x "))
  }
  results_list <- compute_jrqa_measures(rqa_ts_1$RP, rqa_ts_2$RP,
                                        mindiagline = mindiagline,
                                        minvertline = minvertline)
  # Add radius to results_list
  results_list[[1]]$radius <- radius_1_solution
  results_list[[2]]$radius <- radius_2_solution
  results_list
}

# mvjrqa() function  ===========================
#
# This is the main function provided, which will call either
# fixed_radius_mvjrqa() or fixed_recurrence_mvjrqa() depending on
# the value of `setrec`.
#
mvjrqa <- function(ts1,
                   ts2,
                   delay_1,
                   embed_1,
                   radius_1,
                   delay_2,
                   embed_2,
                   radius_2,
                   rescale = 0, # 0 is "none". Allow character argument
                   normalize = 2, # 2 is "z-score". Allow character argument
                   mindiagline = 2,
                   minvertline = 2,
                   tw = 1,
                   whiteline = FALSE,  # TRUE gives error in crqa::tt() internal helper
                   recpt = FALSE,
                   side = "both",
                   metric = "euclidean",
                   datatype = "continuous",
                   setrec = FALSE,
                   targetrec = NULL) {
  # Save how the function was called, so that we can return it. This will
  # facilitate reproducibility.
  command <- match.call()
  # Save the dimensions of the two time series, to be able to return them
  # with the results for improved reproducibility.
  # If the time series is a vector use length(), otherwise use dim()
  get_data_dim <- function(time_series) {
    if (is.null(dim(time_series))) {
      data_dim <- length(time_series)
    } else {
      data_dim <- dim(time_series)
    }
    data_dim
  }
  data_dimensions <- list("ts1" = get_data_dim(ts1),
                          "ts2" = get_data_dim(ts2))
  #
  #  Check whether the radii are fixed or if recurrence rate is fixed
  #
  if (setrec) {
    fixed_rec_result <- fixed_recurrence_mvjrqa(ts1,
                                                ts2,
                                                delay_1,
                                                embed_1,
                                                radius_1,
                                                delay_2,
                                                embed_2,
                                                radius_2,
                                                rescale = rescale,
                                                normalize = normalize,
                                                mindiagline = mindiagline,
                                                minvertline = minvertline,
                                                tw = tw,
                                                whiteline = whiteline,
                                                recpt = recpt,
                                                side = side,
                                                metric = metric,
                                                datatype = datatype,
                                                targetrec = targetrec,
                                                method = method
    )
    # Add how the function was called
    fixed_rec_result$command <- command
    # Add data dimensions
    fixed_rec_result$data_dimensions <- data_dimensions
    return(fixed_rec_result)
  } else {
    fixed_r_result <- fixed_radius_mvjrqa(ts1,
                                          ts2,
                                          delay_1,
                                          embed_1,
                                          radius_1,
                                          delay_2,
                                          embed_2,
                                          radius_2,
                                          rescale = rescale,
                                          normalize = normalize,
                                          mindiagline = mindiagline,
                                          minvertline = minvertline,
                                          tw = tw,
                                          whiteline = whiteline,
                                          recpt = recpt,
                                          side = side,
                                          metric = metric,
                                          datatype = datatype
    )
    # Add how the function was called
    fixed_r_result$command <- command
    # Add data dimensions
    fixed_r_result$datsa_dimensions <- data_dimensions
    return(fixed_r_result)
  }
}
