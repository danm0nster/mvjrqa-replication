## ---------------------------
##
## Script name: 03_empirical_mvjrqa_regression.R
##
## Purpose of script: Perform regression analysis on MvJRQA results
##
## Author: Sebastian Wallot & Dan Moenster
##
## Date Created: August 2024
##
## ---------------------------
##
## Notes:
##
## This script estimates regression models and also produces figures and
## tables related to the models.
##
## ---------------------------

# Wipe work space
rm(list = ls())

# load libraries
library(readr)
library(tidyr)
library(performance)
library(dplyr)
library(lme4)
library(lmerTest)
library(effectsize)
library(modelbased)
library(see)
library(ggplot2)
library(texreg)

source("R/utils.R")
create_dir_if_not_present("tables")

# Perform the regression models for a recurrence rate of 2%
target_rate <- 2
# Accepted maximum deviation from target_RR
rr_tol <- 1.5 

# Load data and remove columns not needed
mvjrqa <- read_csv("results/mvjrqa.csv", show_col_types = FALSE) |>
  select(-c(eye2d_len1, eye2d_len2, eye2d_len3,
            eye3d_len1, eye3d_len2, eye3d_len3)) |> 
  filter(target_RR == target_rate) 

total_obs <- nrow(mvjrqa)
message("Total number of observations: ", total_obs)

# Reshapes data from wide to long:
# splits column names into "EyeType" (2d/3d) and measurement type
mvjrqa_long <- mvjrqa  |>
  pivot_longer(
    cols = c(eye2d_JRR, eye3d_JRR,
             eye2d_DET, eye3d_DET,
             eye2d_LAM, eye3d_LAM,
             eye2d_DENTR, eye3d_DENTR,
             eye2d_RR1, eye3d_RR1,
             eye2d_RR2, eye3d_RR2,
             eye2d_RAD1, eye3d_RAD1,
             eye2d_RAD2, eye3d_RAD2,
             eye2d_file, eye3d_file),
    names_to = c("EyeType", ".value"),
    names_pattern = "eye(\\d+d)_(.*)") |> 
  mutate(
    task_id = factor(task_id),
    subject_id = factor(subject_id),
    EyeType = factor(toupper(EyeType)),
    JRCI = JRR / ((RR1 + RR2) / 2)^2
  )


# Filter on RR1 deviation from target_RR for 2D eye tracking
mvjrqa_1 <- mvjrqa_long |> 
  filter(EyeType == "3D" | (EyeType == "2D" & abs(RR1 - target_RR) <= rr_tol))

# Report number of removed observations from 2D Eye tracking
n_rem_1 <-   nrow(mvjrqa_long) - nrow(mvjrqa_1)
message("Removed ",
        n_rem_1,
        " observations from 2D eye tracking because |RR - RR_target| > ",
        rr_tol)

# Filter on RR1 deviation from target_RR for 3D eye tracking
mvjrqa_clean <- mvjrqa_1 |> 
  filter(EyeType == "2D" | (EyeType == "3D" & abs(RR1 - target_RR) <= rr_tol))

# Report number of removed observations from 3D Eye tracking
n_rem_2 <-  nrow(mvjrqa_1) - nrow(mvjrqa_clean)
message("Removed ",
        n_rem_2,
        " observations from 3D eye tracking because |RR - RR_target| > ",
        rr_tol)

# Estimate the model
m1 <- lmer(JRCI ~ EyeType + (1 | subject_id), data = mvjrqa_clean)

# Uncomment next line to see diagnostic plots
# check_model(m1)
# Print R2 and other model quality measures
model_performance(m1)
r2_nakagawa(m1)
# Standardize the model and summarize it
m1_std <- standardize(m1)
summary(m1_std)
# Just the standardized parameters with 95% CI
standardize_parameters(m1)

# Produce a plot of the model (estimated marginal means) along with the data.
dimension_colours <- c(
  "2D" = "#fd8d3c",
  "3D" = "#b10026"
)

n_2d <- total_obs - n_rem_1
n_3d <- total_obs - n_rem_2

# Compute marginal means for the model
means <- estimate_means(m1, ci = 0.95, by = c("EyeType"))

plot(means) +
  geom_line(data = data.frame(EyeType = c(0.7, 2.3), JRCI = rep(0.01, 2)),
            aes(x = EyeType, y = JRCI),
            colour = "blue") +
  aes(colour = EyeType) +
  geom_point(data = mvjrqa_clean |> filter(target_RR == 2),
             size = 1,
             position = "jitter",
             aes(x = EyeType, y = JRCI, colour = EyeType)) +
  scale_color_manual(values = dimension_colours,
                     name = NULL,
                     aesthetics = c("colour")) +
  ylim(c(0.008, 0.012)) +
  guides(colour = "none",
         colour = guide_legend(position = "top"),
         linestyle = "none") +
  labs(x = "Eye Type",
       y = "JRCI") +
  annotate("text", x = 1, y = 0.012, label = paste("N = ", n_2d), size = 7) +
  annotate("text", x = 2, y = 0.012, label = paste("N = ", n_3d), size = 7) +
  theme_classic() +
  theme(text = element_text(size = 22))


ggsave(filename = "plots/empirical_regression.pdf",
       width = 1.5 * 4, height = 1.5 * 3)

# Tables for main regression and extra models for robustness check

# Robustness check 1: Keep all observations
m2 <- lmer(JRCI ~ EyeType + (1 | subject_id), data =  mvjrqa_long)
model_performance(m2)
# Standardize the model and summarize it
m2_std <- standardize(m2)
summary(m2_std)
# Just the standardized parameters with 95% CI
standardize_parameters(m2)


# Robustness check 2: Remove 3D eye tracking for the trials where 2D is removed
remove_files <- mvjrqa_long |> 
  filter(abs(RR1 - target_RR) > rr_tol) |> 
  pull(`Eye File Name`)

mvjrqa_2 <- mvjrqa_long |> 
  filter(!(`Eye File Name` %in% remove_files))


m3 <- lmer(JRCI ~ EyeType + (1 | subject_id), data =  mvjrqa_2)
model_performance(m3)
# Standardize the model and summarize it
m3_std <- standardize(m3)
summary(m3_std)
# Just the standardized parameters with 95% CI
standardize_parameters(m3)


create_table_models <- function(m_list) {
  # An empty list to hold the models
  table_models <- list()
  # Create four models with different statistics for each model
  # These will become columns in the final table.
  n <- 1
  for (m in m_list) {
    # Model 1 contains the coefficients and SE.
    # Since model 2 also contains SE the idea is to use force.ci = TRUE
    # with model 1 in the call to texreg()
    # Extract results from model
    tr.coef <- texreg::extract(m, include.aic = FALSE, include.bic = FALSE,
                               include.loglik = FALSE)
    # Compute the R-squared of the model
    rsq <-  r2_nakagawa(m)
    # Keep just the first element (Num. obs.) of the gof and add the R-squared
    tr.coef@gof.names <- c(tr.coef@gof.names[1], "Conditional R2")
    tr.coef@gof <- c(tr.coef@gof[1], rsq[["R2_conditional"]])
    tr.coef@gof.decimal <- c(tr.coef@gof.decimal[1], TRUE)
    tr.coef@pvalues <- numeric()
    tr.coef@model.name <- "Estimate"
    table_models[[n]] <- tr.coef
    n <- n + 1
    # Model 2 contains the SE
    tr.se <- texreg::extract(m, include.aic = FALSE, include.bic = FALSE,
                             include.loglik = FALSE)
    tr.se@coef <- tr.se@se
    tr.se@se <- numeric()
    tr.se@pvalues <- numeric()
    tr.se@gof <- numeric()
    tr.se@gof.decimal <- logical()
    tr.se@gof.names <- character()
    tr.se@model.name <- "SE"
    table_models[[n]] <- tr.se
    n <- n + 1
    # Model 3 contains the t-value
    tr.tv <- texreg::extract(m, include.aic = FALSE, include.bic = FALSE,
                             include.loglik = FALSE)
    tr.tv@coef <- tr.coef@coef / tr.tv@se
    tr.tv@se <- numeric()
    tr.tv@pvalues <- numeric()
    tr.tv@model.name <- "$t$"
    tr.tv@gof <- numeric()
    tr.tv@gof.decimal <- logical()
    tr.tv@gof.names <- character()
    table_models[[n]] <- tr.tv
    n <- n + 1
    # Model 4 contains the p-value
    tr.pv <- texreg::extract(m, include.aic = FALSE, include.bic = FALSE,
                             include.loglik = FALSE)
    tr.pv@coef <- tr.pv@pvalues
    tr.pv@se <- numeric()
    tr.pv@pvalues <- numeric()
    tr.pv@model.name <- "$p$"
    tr.pv@gof <- numeric()
    tr.pv@gof.decimal <- logical()
    tr.pv@gof.names <- character()
    table_models[[n]] <- tr.pv
    n <- n + 1
  }
  # Return the list of models
  table_models
}

# Function to create a note, inspired by this thread:
# https://github.com/leifeld/texreg/issues/4
create_model_note <- function(m_list, note_text,
                              width = 0.7,
                              notesize = "small") {
  paste0("\\multicolumn{", length(m_list) + 1,
         "}{p{", width, "\\linewidth}}{\\", notesize, "{", note_text, "}}")
}


m1_models <- create_table_models(list(m1_std))


note_text <- paste(
  "Note: EyeType is a dummy variable that codes for whether 2D or",
  "3D eye movements were used together with the EEG data. Hence, the",
  "predictor EyeType3D provides an estimated increment for the case",
  "where eye movements were 3D with the intercept being the reference",
  "category (i.e., 2D). Values in square brackets are 95\\% CI."
)

m1_note <- create_model_note(m1_models, note_text, width = 0.6)

m1_caption <- paste0(
  "\\textbf{Regression results for MvJRQA.} ", 
  "The table shows fixed effects as standardized regression coefficients ",
  "for the regression with JRCI as the dependent variable computed for EEG ",
  "coupled with either 2D or 3D eye movement dynamics. ",
  "Subsystem recurrence was fixed at $RR = ", target_rate, "\\%$ and ",
  n_rem_1, " observations were removed because the recurrence rate differed ",
  "from the target rate by more than ", rr_tol, " percentage points."
)

# Print a simple version to screen
texreg::screenreg(m1_models,
                  digits = 4,
                  stars = 0, # Omit stars for p-values
                  ci.test = NA, # Omit stars for CI test
                  ci.force = c(TRUE, FALSE, FALSE, FALSE),
                  custom.note = ""
)

# This one line using sjPlots, basically does the same thing,
# but does not save to LaTeX and has CI on same line (too wide)
# Install sjPlots and uncomment to use
# sjPlot::tab_model(m1_std, show.se = TRUE, show.stat = TRUE)

# Save the full version to file in LaTeX format
texreg::texreg(m1_models,
               stars = 0, # Omit stars for p-values
               ci.test = NA, # Omit stars for CI test
               ci.force = c(TRUE, FALSE, FALSE, FALSE),
               caption.above = TRUE,
               caption = m1_caption,
               label = "tab:regression",
               custom.note = m1_note,
               file = "tables/regression_main.tex"
)

# Extra models for robustness checks

robutness_caption <- paste(
  "\\textbf{Robustness check on regression results for MvJRQA.}",
  "Model 1 is the same model shown in Table~\\ref{tab:regression},",
  "but including all", 2 * total_obs, "observations, i.e., no,",
  "observations were removed, even if $RR - RR_T > \\delta$. Model 2",
  "represents the other extreme, viz.\\ removing the MvJRQA results",
  "for EEG combined with both the 2D and the 3D eye tracking results",
  "in all the cases where $RR - RR_T > \\delta$ for 2D eye tracking."
)

# Prepare the standardized models for a texreg table
std_models <- list(m2_std, m3_std)
all_models <- create_table_models(std_models)
texreg::screenreg(all_models,
                  custom.header = list("Model 1" = 1:4, "Model 2" = 5:8),
                  stars = 0, # Omit stars for p-values
                  single.row = FALSE, # CI after coef on same row id TRUE
                  ci.force = rep(c(TRUE, FALSE, FALSE, FALSE),
                                 length(std_models)),
                  ci.test = NA, # Do not add stars for CI test
                  custom.note = ""
)

# sjPlots version.
# Install sjPlots and uncomment to use
# sjPlot::tab_model(m2_std, m3_std, show.se = TRUE, show.stat = TRUE)

texreg::texreg(all_models,
               custom.header = list("Model 1" = 1:4, "Model 2" = 5:8),
               stars = 0, # Omit stars for p-values
               single.row = FALSE, # CI after coef on same row id TRUE
               ci.force = rep(c(TRUE, FALSE, FALSE, FALSE),
                              length(std_models)),
               ci.test = NA, # Do not add stars for CI test
               custom.note = "",
               caption.above = TRUE,
               caption = robutness_caption,
               label = "tab:robustness",
               file = "tables/regression_robustness.tex"
)

# Report mean and SD for the RR used for the main analysis (m1)
rr_vals <- c(mvjrqa_clean$RR1, mvjrqa_clean$RR2)
rr_mean <- mean(rr_vals)
rr_sd <- sd(rr_vals)
message("RR M =  ", round(rr_mean, 2))
message("RR SD = ", round(rr_sd, 2))
