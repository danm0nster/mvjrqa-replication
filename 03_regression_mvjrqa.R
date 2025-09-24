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
library(car)
library(dplyr)
library(lme4)
library(lmerTest)

# Load data and remove columns not needed
df <- read_csv("results/mvjrqa.csv", show_col_types = FALSE) |>
  select(-c(eye2d_len1, eye2d_len2, eye2d_len3,
            eye3d_len1, eye3d_len2, eye3d_len3))

# Reshapes data from wide to long:
# splits column names into "Eye_Type" (2d/3d) and measurement type
df2 <- df  |>
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
    names_pattern = "eye(\\d+d)_(.*)"
  )

# factor and log-transform variables
df2 <- df2 |>
  mutate(
    task_id = factor(task_id),
    subject_id = factor(subject_id),
    Eye_Type = factor(Eye_Type),
    logJRR = log(JRR + 0.1),
    logDET = log(DET + 0.1),
    logDENTR = log(DENTR + 0.1)
  )

# Models of eye type on log of JRR, DET, and DENTR
m1 <- lmer(logJRR ~ Eye_Type + (1 | subject_id), data = df2)
summary(m1)
vif(m1)
r2_nakagawa(m1)

m1 <- lmer(logDET ~ Eye_Type + (1  | subject_id), data = df2)
summary(m1)
vif(m1)
r2_nakagawa(m1)

m1 <- lmer(logDENTR ~ Eye_Type + (1 | subject_id), data = df2)
summary(m1)
vif(m1)
r2_nakagawa(m1)


# Do analysis of JRCI
df2$JRCI <-  df2$JRR / (rowMeans(cbind(df2$RR1, df2$RR2))^2)
m1 <- lmer(JRCI ~ Eye_Type + (1 | subject_id), data = df2)
summary(m1)
vif(m1)
r2_nakagawa(m1)
