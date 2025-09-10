## ---------------------------
##
## Script name: 02_merge_mvjrqa_files.R
##
## Purpose of script: Gather and merge files from MvJRQA analysis
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

# wipe workspace
rm(list=ls())

# load libraries
library(dplyr)
library(readr)

### LOAD AND CONCATENATE RESUTLS
eye2d_path <- "results/eye2d"
eye3d_path <- "results/eye3d"

files = list.files(path = eye2d_path, pattern = "results\\.txt$", full.names = TRUE)
eye2d = data.frame()
for (f in files) {
  df = read.table(f, header = TRUE, sep = ",")
  df$eye2d_file = f
  eye2d = rbind(eye2d, df)
}
# Save the results to a file
write.table(eye2d, file = "results/eye2d.txt", sep = ",", row.names = FALSE)

files = list.files(path = eye3d_path, pattern = "results\\.txt$", full.names = TRUE)
eye3d = data.frame()
for (f in files) {
  df = read.table(f, header = TRUE, sep = ",")
  df$eye3d_file = f
  eye3d = rbind(eye3d, df)
}
write.table(eye3d, file = "results/eye3d.txt", sep = ",", row.names = FALSE)

# Exclude these two files, since there are issues with them
eeg_exclude <- c("21_3_4.edf", "21_3_5.edf")

# Read performance data and rename some variables for joining
perf <- read_csv(
  "empirical_data/PerformanceScores.csv",
  quote = "'",
  show_col_types = FALSE) |> 
  rename("Gender(F:Female, M:Male)" = "Gender(F:Female; M:Male)") |> 
  filter(!(`EEG File Name` %in% eeg_exclude)) |> 
  rename(
    subject_id = `Subject ID`,
    task_id = `Task ID`,
    try = Try
  )

# eye2d <- read_csv("results/eye2d.txt")
# eye3d <- read_csv("results/eye3d.txt")

# Merge performance data with the MvJRQA results
merged_data <- left_join(perf, eye2d, by = c("subject_id", "task_id", "try"))
merged_data <- left_join(merged_data, eye3d, by = c("subject_id", "task_id", "try"))

# Write result to CSV file
# TODO: Use write_csv()
write_csv(merged_data, file = "results/mvjrqa.csv")
