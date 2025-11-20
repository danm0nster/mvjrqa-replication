## ---------------------------
##
## Script name: 03_regression_mvjrqa.R
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
##
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

# Perform the regression models for a recurrence rate of 2%
target_rate <- 2

# Load data and remove columns not needed
mvjrqa <- read_csv("results/mvjrqa.csv", show_col_types = FALSE) |>
  select(-c(eye2d_len1, eye2d_len2, eye2d_len3,
            eye3d_len1, eye3d_len2, eye3d_len3)) |> 
  filter(target_RR == target_rate) 

# Reshapes data from wide to long:
# splits column names into "Eye_Type" (2d/3d) and measurement type
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
    names_to = c("Eye_Type", ".value"),
    names_pattern = "eye(\\d+d)_(.*)")
    
total_obs <- nrow(mvjrqa_long)

rr_tol <- 1.5 # Accepted maximum deviation from target_RR

# Filter on RR1 deviation from target_RR
mvjrqa_1 <- mvjrqa_long |> 
  filter(abs(RR1 - target_RR) <= rr_tol)

# Report number of removed observations from 2D Eye tracking
n_rem_1 <-  total_obs - nrow(mvjrqa_1)
message("Observations removed because RR1 is outside tolerance: ", n_rem_1)

# Filter on RR1 deviation from target_RR
mvjrqa_2 <- mvjrqa_1 |> 
  filter(abs(RR2 - target_RR) <= rr_tol)

# Report number of removed observations from 3D Eye tracking
n_rem_2 <-  nrow(mvjrqa_1) - nrow(mvjrqa_2)
message("Observations removed because RR2 is outside tolerance: ", n_rem_2)

# TODO: make a record of how many observations we remove because
# we require RR1 and RR2 not to deviate too much from target_RR.

# factor and log-transform variables
mvjrqa_clean <- mvjrqa_2 |>
  mutate(
    task_id = factor(task_id),
    subject_id = factor(subject_id),
    Eye_Type = factor(Eye_Type),
    logJRR = log(JRR + 0.1),
    logDET = log(DET + 0.1),
    logDENTR = log(DENTR + 0.1)
  )

# Models of eye type on log of JRR, DET, and DENTR
m1 <- lmer(logJRR ~ Eye_Type + (1 | subject_id), data = mvjrqa_clean)
# summary(m1)
r2_nakagawa(m1)

m1_std <- standardize(m1)
summary(m1_std)

m2 <- lmer(logDET ~ Eye_Type + (1  | subject_id), data = mvjrqa_clean)
# summary(m2)
r2_nakagawa(m2)

m2_std <- standardize(m2)
summary(m2_std)

m3 <- lmer(logDENTR ~ Eye_Type + (1 | subject_id), data = mvjrqa_clean)
# summary(m3)
r2_nakagawa(m3)

m3_std <- standardize(m3)
summary(m3_std)


# Do analysis of JRCI
mvjrqa_clean$JRCI <-  mvjrqa_clean$JRR /
  (rowMeans(cbind(mvjrqa_clean$RR1, mvjrqa_clean$RR2))^2)
m4 <- lmer(JRCI ~ Eye_Type + (1 | subject_id), data = mvjrqa_clean)
# summary(m4)
r2_nakagawa(m4)

m4_std <- standardize(m4)
summary(m4_std)


