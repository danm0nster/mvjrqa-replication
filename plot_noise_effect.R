## ---------------------------
##
## Script name: plot_noise_effect.R
##
## Purpose of script: Make a plot to illustrate the effect of noise
##
## Author: Dan Moenster & Sebastian Wallot
##
## Date Created: May 6 2025
##
## ---------------------------
##
## Notes:
## 
## The script only shows the effect of observational noise, not dynamical noise.
##
## The Lorenz system with the canonical parameters is used to show the effect
## of adding noise to a signal. Random initial conditions are used.
##
## The column `coupling` is used to mean noise.
##
## ---------------------------

source("R/lorenz_harmonic.R")
source("R/utils.R")
source("R/plot_utils.R")
source("R/mvjrqa.R")

set.seed(29548)
x_init <- runif(1, min = -10, max = 10)
y_init <- runif(1, min = -10, max = 10)
z_init <- runif(1, min = 0, max = 40)
# Here u and v are not used, but the model computes them, so provide
# initial values.
u_init <- runif(1, min = -10, max = 10)
v_init <- runif(1, min = -10, max = 10)

sample_length <- 500
sample_skip <- 100

model_1 <- lorenz_harmonic(n = sample_length,
                              skip = sample_skip,
                              coupling = 0, # We do not need coupled systems
                              initial_state = c(x = x_init,
                                                y = y_init,
                                                z = z_init,
                                                u = u_init,
                                                v = v_init))


ts_1 <- extract_lorenz(model_1)


# Now perform mvjrqa for different fixed RR values
rr_values <- c(0.2, 0.3, 0.6, 1.25, 2.5, 5, 10, 20, 40, 80)

compute_jrc <- function(ts_1, ts_2,
                        target_rr = NULL,
                        delay = 1,
                        embed = 1) {
  joint_rr <- data.frame()
  for (rr_target in target_rr) {
    mvjrqa_result <- mvjrqa(ts_1,
                            ts_2, 
                            delay_1 = delay,
                            embed_1 = embed,
                            radius_1 = NA,
                            delay_2 = delay,
                            embed_2 = embed,
                            radius_2 = NA,
                            setrec = TRUE,
                            targetrec = rr_target)
    rqa_system_1 <- as.data.frame(mvjrqa_result[[1]]) |> 
      mutate(system = "1", 
             coupling = 0, # Set coupling to zero
             RR_target = rr_target)
    rqa_system_2 <- as.data.frame(mvjrqa_result[[2]]) |> 
      mutate(system = "2", 
             coupling = 0, # Set coupling to zero
             RR_target = rr_target)
    rqa_joint <- as.data.frame(mvjrqa_result[[3]]) |> 
      mutate(system = "j", 
             coupling = 0, # Set coupling to zero
             RR_target = rr_target)
    # Return the results in one data frame
    joint_rr <- bind_rows(joint_rr,
                          rqa_system_1,
                          rqa_system_2,
                          rqa_joint)
  }
  return(joint_rr)
}

# Function to add noise
add_noise <- function(data_matrix, noise_level = 0.05) {
  noisy_matrix <- data_matrix
  # Add Gaussian noise with sd = noise_level * sd(signal) for each column
  for (j in 1:ncol(data_matrix)) {
    signal_sd <- sd(data_matrix[, j])
    noise <- rnorm(n = nrow(data_matrix),
                   mean = 0,
                   sd = noise_level * signal_sd)
    noisy_matrix[, j] = noisy_matrix[, j] + noise
  }
  return(noisy_matrix)
}

#
# Look at identical signals with varying levels of noise
#

noise_levels <- c(0, 0.01, 0.05, 0.10, 1)

noise_data <- data.frame()
for (nl in noise_levels) {
  # Add noise to two versions of ts_1.
  ts_1_noise_nl_a <- add_noise(ts_1, noise_level = nl)
  ts_1_noise_nl_b <- add_noise(ts_1, noise_level = nl)
  ts_1_nl_result <-  compute_jrc(ts_1_noise_nl_a, ts_1_noise_nl_b,
                                 target_rr = rr_values)
  # Note: This is a slight abuse of the coupling variable to mean noise.
  ts_1_nl_result$coupling <- nl
  noise_data <- bind_rows(noise_data, ts_1_nl_result)
}

jrc_plot(noise_data, errorbars = "Off", coupling_name = "xi",
         linewidth = 0.8, linetype = "dashed") +
  scale_color_manual(values = rep("black", length(noise_levels)),
                     aesthetics = c("colour", "fill")) 

plot_file_name <- paste0("plots/lorenz_noise.pdf")
ggsave(plot_file_name, width = 1.5 * 4, height = 1.5 * 3)
