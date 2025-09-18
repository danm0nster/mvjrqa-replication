## ---------------------------
##
## File name: lorenz96.R
##
## Purpose:  Numerically solve Lorenz 96 model coupled to harmonic oscillator.
##
## Author: Dan Moenster & Sebastian Wallot
##
## Date Created: May 27 2025
##
## ---------------------------
##
## Notes:
##
## The Lorenz 96 model can have a variable number of dimensions (variables).
## This can be set using the n_var argument which defaults to 5.
##
## Also contains helper functions to extract the results for each of the
## two sub-systems:
##  extract_lorenz96
##  extract_oscillator
##
## ---------------------------

library(deSolve)
library(data.table)
library(dplyr)

#
# Lorenz 96 model coupled to harmonic oscillator.
#
lorenz96_harmonic <- function(n = 1000,
                              skip = 100,
                              dt = 0.01,
                              n_var = 5,
                              forcing = 8,
                              coupling = 0.5,
                              initial_state = list(x = c(1.01, rep(1, 4)),
                                                   u = 10,
                                                   v = 0)) {
  # Check that the initial state of the Lorenz 96 system has the same number
  # of values as the number of variables. Otherwise issue a warning and set
  # the initial state to (1, 1, 1, 1, ...) with a small random perturbation.
  # NOTE: Using lengths() not length() to get lengths of list elements
  x_length <- lengths(initial_state["x"])
  if (x_length != n_var) {
    warning(
      paste0("Number of x-values in inital state (",
            x_length,
            ") differs from n_var (",
            n_var,
            ").\n",
            "    Initial values will be random.")
    )
    initial_state$x <- rnorm(n_var, mean = 1, sd = 0.01)
  }
  # Internal function to compute derivatives for the ode solver.
  # x: The Lorenz 96 system
  # u, v: position and velocity of harmonic oscillator, perturbed by Lorenz-x[1]
  # A force proportional to square of the x[1]-component of the Lorenz 96 model
  # is exerted on the oscillator with a coupling constant chi.
  coupled <- function(t, state, parameters) {
    x_indices <- grep("x", labels(state))
    x <- state[x_indices]
    with(as.list(c(state, parameters)), {
      x_shifted_pos_1 <- data.table::shift(x, -1, type = "cyclic")
      x_shifted_neg_1 <- data.table::shift(x, +1, type = "cyclic")
      x_shifted_neg_2 <- data.table::shift(x, +2, type = "cyclic")
      dx <- (x_shifted_pos_1 - x_shifted_neg_2) * x_shifted_neg_1 - x + f
      du <- v
      dv <- -k * u + c * x[1]**2
      return(list(c(dx, du, dv)))
    })
  }
  # Parameters for the model. The constant c is taken from the function
  # call to lorenz_harmonic, so it can be varied.
  params <- c(f = forcing,
              k = 10,
              c = coupling)
  # A vector of time values created from the step-size, dt, the
  # number of points to skip (to avoid transient behavior) and
  # the number of points to output, n.
  times <- seq(0, dt * (skip + n), dt)
  # Use ode to solve the system
  ode_data <- ode(y = unlist(initial_state), times = times,
                  func = coupled, parms = params)
  # Convert the data to a data frame
  model_data <- as.data.frame(ode_data)
  if (skip > 0) {
    model_data <- model_data[(skip + 1):(skip + n), ]
  }
  return(model_data)
}

# Function to extract the Lorenz part of the data as a matrix
# The input model_df is the output from lorenz96_harmonic()
extract_lorenz96 <- function(model_df) {
  return(
    as.matrix(
      model_df |> select(starts_with("x"))
    )
  )
}

# Function to extract the oscillator part of the data as a matrix
# The input model_df is the output from lorenz96_harmonic()
extract_oscillator <- function(model_df) {
  return(
    as.matrix(
      model_df |> select(u, v)
    )
  )
}
