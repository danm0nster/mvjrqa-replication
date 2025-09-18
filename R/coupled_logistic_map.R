## ---------------------------
##
## File name: coupled_logistic_map.R
##
## Purpose:  Numerically solve coupled logistic maps driven by cosine.
##
## Author: Dan Moenster & Sebastian Wallot
##
## Date Created: Apr 12 2025
##
## ---------------------------
##
## Notes:
##
## Since this is a discrete time model it is "solved" by iteration.
##
## Also contains helper functions to extract the results for each of the
## two sub-systems (X and Y) as well as the driver (H):
##  extract_X
##  extract_Y
##  extract_H
##
## ---------------------------
library(dplyr)

coupled_logistic_map <- function(x0 = 0.2,
                                 y0 = 0.6,
                                 rx = 3.65,
                                 ry = 3.8,
                                 bxy = 0,
                                 byx = 0.4,
                                 etax = 0,
                                 etay = 0,
                                 p = 20,
                                 phi = 0,
                                 N = 100,
                                 N_skip = 0) {
  if (p <= 0) {
    stop("The period, p, must be positive. Has value: ", p)
  }
  X <- rep(0, N)
  Y <- rep(0, N)
  H <- rep(0, N)
  # Time that ticks during the transient N_skip iterations is needed for H(t),
  # since the function needs "absolute" time.
  time <- 0
  if (N_skip > 0) {
    #
    # Iterate for N_trans generations without collecting data (only
    # X[1] and Y[1] are updated, so the last data point will be
    # used as the new X[1] and Y[1]).
    #
    H0 <- cos(phi) # time = 0, so cos(2 * pi * time / p + phi) = cos(phi)
    X0 <- x0
    Y0 <- y0
    for (t in 1:N_skip) {
      time <- time + 1
      H[1] <- cos(2 * pi * time / p + phi)
      X[1] <- X0 * ((rx + etax * H0) * (1 - X0) - bxy * Y0)
      Y[1] <- Y0 * ((ry + etay * H0) * (1 - Y0) - byx * X0)
      H0 <- H[1]
      X0 <- X[1]
      Y0 <- Y[1]
    }
  } else {
    H0 <- 0
    X0 <- x0
    Y0 <- y0
    time <- time + 1
    H[1] <- cos(2 * pi * time / p + phi)
    X[1] <- X0 * ((rx + etax * H0) * (1 - X0) - bxy * Y0)
    Y[1] <- Y0 * ((ry + etay * H0) * (1 - Y0) - byx * X0)
  }

  # Iterate the coupled maps for N generations
  for (t in 2:N) {
    time <- time + 1
    H[t] <- cos(2 * pi * time / p + phi)
    X[t] <- X[t - 1] * ((rx + etax * H[t - 1]) * (1 - X[t - 1]) -
                          bxy * Y[t - 1])
    Y[t] <- Y[t - 1] * ((ry + etay * H[t - 1]) * (1 - Y[t - 1]) -
                          byx * X[t - 1])
  }

  return(
    data.frame(
      time = 1:N,
      X = X,
      Y = Y,
      H = H
    )
  )
}


# Functions to extract, respectively, the X, Y and H variables of the data
# as a vector. The input model_df is the output from coupled_logistic_map()

extract_X <- function(model_df) {
  return(model_df |> pull(X))
}

extract_Y <- function(model_df) {
  return(model_df |> pull(Y))
}

extract_H <- function(model_df) {
  return(model_df |> pull(H))
}
