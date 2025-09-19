## ---------------------------
##
## Script name: create_mvjrp_animation.R
##
## Purpose of script: Create animation movie to illustrate MvJRQA.
##
## Author: Dan Moenster & Sebastian Wallot
##
## Date Created: Aug 29 2025
##
## ---------------------------
##
## Notes:
##
## This script will produce an animated movie showing an example of how joint
## recurrence rate and the recurrence plot of the driven system change when
## the coupling between the two sub-systems is changed.
##
## ---------------------------

library(ggplot2)
library(cowplot)
library(purrr)
library(animation)
library(gridExtra)
library(crqa)

source("R/mvjrqa.R")
source("R/lorenz_harmonic.R")
source("R/utils.R")

# Make sure directory "./plots" exists. This is where movie and plot will be
# saved.
create_dir_if_not_present("plots")

# None of the time series need to be further embedded,
# so delay and embed are set to 1.
delay <- 1
embed <- 1
radius <- 0
fixed_rr <- 2

# Initialise lists and data frames to hold the data
rp_1_list <- list()
rp_2_list <- list()
rp_j_list <- list()
rp_1_matrix <- list()
rp_2_matrix <- list()
rp_j_matrix <- list()
lagged_rr_plot <- list()
time_series_data <- data.frame()
joint_rr <- data.frame()
lagged_rr_data <- data.frame()

# Loop over values of the coupling constant, cc.
# couplings <- seq(0, 0.38, 0.02)
couplings <- seq(0, 1.28, 0.04)


# TODO: Put data generation into a separate script
for (cc in couplings) {
  # Run the model to produce data
  model_data <- lorenz_harmonic(n = 500,
                                skip = 500,
                                coupling = cc,
                                initial_state = c(x = 10,
                                                  y = 10,
                                                  z = 10,
                                                  u = 10,
                                                  v = 0)
                                )
  
  # Add coupling to model_data and save
  time_series_data <- bind_rows(
    time_series_data,
    model_data |> mutate(c = cc)
  )
  
  # Extract the two time series from model_data
  ts_1 <- extract_lorenz(model_data)
  ts_2 <- extract_oscillator(model_data)
  
  rqa_results <- mvjrqa(ts_1, ts_2,
                        delay, embed, radius,
                        delay, embed, radius,
                        setrec = TRUE,
                        targetrec = fixed_rr)
  # TODO: Remove this test code
  # Test
  # plot_rp(rqa_results$rp_list[[1]])

  # TODO: Is this needed?
  # Save the sparse matrix versions
  rp_1_matrix[[length(rp_1_matrix) + 1]] <- rqa_results$rp_list[[1]]
  rp_2_matrix[[length(rp_2_matrix) + 1]] <- rqa_results$rp_list[[2]]
  rp_j_matrix[[length(rp_j_matrix) + 1]] <- rqa_results$rp_list[[3]]
  
  # Save the RR in a data frame for easy access
  joint_rr <- rbind(joint_rr,
                    data.frame(c = cc,
                               rr_1 = rqa_results[[1]]$RR,
                               rr_2 = rqa_results[[2]]$RR,
                               rr_j = rqa_results[[3]]$RR)
                    )
  
  rp_1_list[[length(rp_1_list) + 1]] <- plot_rp(rqa_results$rp_list[[1]])
  rp_2_list[[length(rp_2_list) + 1]] <- plot_rp(rqa_results$rp_list[[2]])
  rp_j_list[[length(rp_j_list) + 1]] <- plot_rp(rqa_results$rp_list[[3]])
} # End loop over coupling


# Note this zooms in on the first 100 x 100 points in the RP.
# TODO: Set a parameter to change this. Maybe no zoom at all?
plot_size <- dim(rp_1_matrix[[1]])[1] # Can be set to smaller value to zoom in

###########################
# Look at part of one JRP #
###########################
animation_ps2 <- list()
animation_plots_jrp <- list()
animation_plots_rp2 <- list()
animation_plots_all <- list()

ani.options(
  ani.width = 720,
  ani.height = 720,
  auto.browse = FALSE,
  autoplay = FALSE
)

for (p in 1:length(rp_j_list)) {
  # Use cowplot::ggdraw() + theme() to set background to white.
  
  animation_ps2[[p]] <-  cowplot::ggdraw(
    ggplot(data = time_series_data |> filter(c == joint_rr$c[p]),
           aes(x = u, y = v)) +
      geom_path() +
      theme_classic() +
      xlim(c(min(time_series_data$u), max(time_series_data$u))) + 
      ylim(c(min(time_series_data$v), max(time_series_data$v))) + 
      annotate("text", 
               x = min(time_series_data$u) + 2,
               y = max(time_series_data$v), 
               hjust = "left", 
               label = paste("c = ", round(joint_rr$c[p], 2))) +
      labs(title = "Oscillator trajectory") +
      theme(aspect.ratio = 1)
  ) +
    theme(plot.background = element_rect(fill="white", color = NA))
  
  animation_plots_rp2[[p]] <-  cowplot::ggdraw(
    rp_2_list[[p]] + 
      xlim(c(1, plot_size)) + 
      ylim(c(1, plot_size)) +
      annotate("label", 
               x = 10, 
               y = 0.95 * plot_size,
               fill = "white",
               label.size = 0,
               hjust = "left",
               label = paste("c = ", round(joint_rr$c[p], 2),
                             "\nRR = ", round(joint_rr$rr_2[p], 2))) +
      labs(title = "Oscillator RP")
  ) +
    theme(plot.background = element_rect(fill="white", color = NA))
  animation_plots_jrp[[p]] <-  cowplot::ggdraw(
    rp_j_list[[p]] + 
      xlim(c(1, plot_size)) + 
      ylim(c(1, plot_size)) +
      annotate("text",
               x = 10,
               y = 0.95 * plot_size,
               hjust = "left",
               label = paste("c = ", round(joint_rr$c[p], 2),
                             "\nRR = ", round(joint_rr$rr_j[p], 2))) +
      labs(title = "Joint RP")
  ) +
    theme(plot.background = element_rect(fill="white", color = NA))
  
  animation_plots_all[[p]] <- arrangeGrob(
    grobs = list(
      animation_ps2[[p]],
      animation_plots_rp2[[p]],
      animation_plots_jrp[[p]],
      ggplot(joint_rr |> filter(c <= joint_rr$c[p]),
             aes(x = c, y = rr_j)) +
        geom_point() +
        geom_line() +
        xlab("coupling strength") +
        ylab("Joint RR") +
        expand_limits(x = max(joint_rr$c),
                      y = max(joint_rr$rr_j)) +
        theme_classic()
    ),
    # heights = unit(1, "npc"),
    # widths = unit(c(0.5, 0.5), "npc"),
    layout_matrix = matrix(c(1, 2, 3, 4, 4, 4), ncol = 3, byrow = TRUE),
    nrow = 2, ncol = 3
  )
}


animation::saveVideo(
  expr = {
    purrr::walk(
      animation_plots_all,
      ~ plot(.)
    )
  },
  video.name = "plots/mvjrp_animation.mp4", # Movie_S1 for the paper
  other.opts = "-s 1080x1080 -pix_fmt yuv420p -preset veryslow -crf 0"
)

# Save final frame in PDF
ggsave(
  filename = "plots/mvjrp_animation_snapshot.pdf",
  plot = animation_plots_all[[length(animation_plots_all)]],
  width = 8,
  height = 8
)

