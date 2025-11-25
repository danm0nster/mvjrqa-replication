## ---------------------------
##
## Script name: 01_empirical_mvjrqa_analysis.R
##
## Purpose of script: Perform MvJRQA analysis in parallel
##
## Author: Sebastian Wallot & Dan Moenster
##
## Date Created: August 2024
##
## ---------------------------
##
## Notes:
##
## This script uses the parallel package for parallel processing.
## Please set the number of worker processes you wish to use by changing
## the value of the variable num_workers below. It should be less than the
## number of CPU cores available on your system.
##
## Note also that around 125 GB of RAM per core is needed, so you need
## sufficient memory for the number of cores used. E.g., if you have 256 GB
## of RAM, you can only use two cores (workers).
##
## In order to run this script, you must first download the empirical data set
## by Shafei & Shadpour (2023) from Physionet.
##
## Link: https://doi.org/10.13026/kyjw-p786
##
## Note that this is a large data set with a size of 7.8 GB.
##
## After downloading, place the contents in the directory called
## "empirical_data" that is a sub-directory of the directory where this
## R script is located.
##
## Note that there is already a file in the directory called
## PerformanceScores_fixed.csv, where we have fixed some errors that were
## present in the original file at the time of publication.
##
## Schematically is should look like below, where "..." means some files
## have been left out.
##
## .
## ├── 01_mvjrqa_analysis_empirical.R
## │   ...
## ├── empirical_data
## │   ├── EEG_FLS
## │   │   ├── 10_1_1.edf
## │   │   ├── 10_1_2.edf
## │   │   ...
## │   ├── EYE_FLS
## │   │   ├── 10_1_1.csv
## │   │   ├── 10_1_2.csv
## │   │   ...
## │   ├── LICENSE.txt
## │   ├── PerformanceScores.csv
## │   ├── PerformanceScores_fixed.csv
## │   ├── RECORDS
##
## ---------------------------

# wipe workspace
rm(list = ls())

# load libraries and source functions
library(readr)
library(dplyr)
library(edf)
source("R/mvjrqa.R")
source("R/utils.R")
library(parallel)

# Set the number of workers (# worker processes <= # CPU cores)
num_workers <- 1

# Individual files are written to results/eye2d and results/eye3d
eye2d_path <- "results/eye2d"
eye3d_path <- "results/eye3d"

# Create directories folders "results/eye2d" and "results/eye3d"
# if they do not exist
create_dir_if_not_present(eye2d_path)
create_dir_if_not_present(eye3d_path)

# Exclude the two EEG files where the two eye files that are not present.
eeg_exclude <- c("21_3_4.edf", "21_3_5.edf")
# Exclude some large files, to avoid integer overflow in vectors and matrices.
large_files <- c("24_3_4.edf", "12_3_1.edf", "12_3_2.edf", "3_1_1.edf",
                 "20_3_1.edf", "20_3_2.edf", "16_3_1.edf")
eeg_exclude <- c(eeg_exclude, large_files)

# Get the file names from the performance data file
# NOTE: PerformanceScores.csv from PhysioNet contains errors.
# We therefore use a fixed version.
perf_file <- "empirical_data/PerformanceScores_fixed.csv"
perf <- read_csv(perf_file,
                 # Only for the PhysioNet version, not the fixed version
                 # quote = "'", 
                 show_col_types = FALSE) |>
  filter(!(`EEG File Name` %in% eeg_exclude))

eeg_data <- perf$`EEG File Name`
eye_data <- perf$`Eye File Name`

data_files <- data.frame(
  eeg_data = perf$`EEG File Name`,
  eye_data = perf$`Eye File Name`,
  subject_id = perf$`Subject ID`,
  task_id = perf$`Task ID`,
  try = perf$Try
)

# Create a vector of row indices for parallelization
i_vec <-  1:nrow(data_files) # vector defining the data files for loop_data_fun

### Function for parallelization
loop_data_fun <-  function(i, x) {
  targetRR <- c(1, 2, 4, 8, 16, 32, 64) # The RR values to compute
  # load and sort eeg data
  eeg <- read.edf(paste0("empirical_data/EEG_FLS/", data_files$eeg_data[i]),
                  read.annotations = TRUE, header.only = FALSE)
  eeg <- data.frame(eeg$signal) # convert to data frame
  # remove t-variables and downsample to 50 Hz
  eeg <- eeg[seq(1, dim(eeg)[1], by = 10), seq(1, dim(eeg)[2], by = 2)]
  drops <- c("EEGHEOGRCPz.data",
             "EEGHEOGLCPz.data",
             "EEGVEOGUCPz.data",
             "EEGVEOGLCPz.data") # list EOGs to drop from data frame
  eeg <- eeg[, !(names(eeg) %in% drops)] # drop EOGs from data frame
  rm(list = "drops")
  # load eye data
  df_eye <- read.csv(paste0("empirical_data/EYE_FLS/", data_files$eye_data[i]))
  eye_2d <- data.frame(df_eye[, 1:2]) # extract 2d eye gaze
  eye_3d <- data.frame(df_eye[, 3:5]) # extract 3d gaze
  # cut off superfluous data points
  if (dim(eeg)[1] > dim(df_eye)[1]) {
    eeg <- eeg[1:dim(eye_2d)[1], ] # cut eeg record
  } else if (dim(eeg)[1] < dim(df_eye)[1]) {
    eye_2d <- eye_2d[1:dim(eeg)[1], ]# cut eye records
    eye_3d <- eye_3d[1:dim(eeg)[1], ]# cut eye records
  } else {
    # do nothing
  }
  le1 <- dim(eeg)[1]
  rm(list = "df_eye")
  # remove NAs from all data sets
  naCount <- rowSums(is.na(eeg)) +
    rowSums(is.na(eye_2d)) +
    rowSums(is.na(eye_3d))
  eeg <- eeg[naCount == 0, ]
  eye_2d <- eye_2d[naCount == 0, ]
  eye_3d <- eye_3d[naCount == 0, ]
  rm(list = "naCount")
  le2 <- dim(eeg)[1]
  # Calculate difference values for each data set
  eeg <- data.frame(lapply(eeg, diff))
  eye_2d <- data.frame(lapply(eye_2d, diff))
  eye_3d <- data.frame(lapply(eye_3d, diff))
  # Calculate means and SDs for EEG data
  stats <- sapply(eeg, function(x) c(mean = mean(x), sd = sd(x)))
  # remove variables with 0 SD
  eeg <- eeg[, stats[2, ] > 0]
  stats <- stats[, stats[2, ] > 0]
  rm(list = c("stats"))
  # length of final data set
  le3 <- dim(eeg)[1]
  # run analysis on eeg and 2d eye data
  eeg_eye2d <- data.frame()
  eeg_eye3d <- data.frame()
  for (rr_target in targetRR) {
    res <- mvjrqa(as.matrix(eye_2d),
                  as.matrix(eeg),
                  delay_1 = 1,
                  embed_1 = 1,
                  radius_1 = 10,
                  delay_2 = 1,
                  embed_2 = 1,
                  radius_2 = 1,
                  rescale = 0,
                  normalize = 0,
                  mindiagline = 2,
                  minvertline = 2,
                  tw = 1,
                  whiteline = TRUE,
                  recpt = FALSE,
                  side = "both",
                  metric = "euclidean",
                  datatype = "continuous",
                  setrec = TRUE,
                  targetrec = rr_target)
    # store RR for individual RPs and JRP
    eeg_eye2d <- rbind(
      eeg_eye2d,
      data.frame(
        subject_id = data_files$subject_id[i],
        task_id = data_files$task_id[i],
        try = data_files$try[i],
        target_RR = rr_target,
        eye2d_RR1 = res[[1]]$RR,
        eye2d_RR2 = res[[2]]$RR,
        eye2d_JRR = res[[3]]$RR,
        eye2d_DET = res[[3]]$DET,
        eye2d_LAM = res[[3]]$LAM,
        eye2d_DENTR = res[[3]]$DENTR,
        eye2d_RAD1 = res[[1]]$radius,
        eye2d_RAD2 = res[[2]]$radius,
        eye2d_len1 = le1,
        eye2d_len2 = le2,
        eye2d_len3 = le3
      )
    )
  write.table(eeg_eye2d,
              file = paste0(eye2d_path, "/eye2d_", i, "_results.txt"),
              sep = ",", row.names = FALSE)
  rm(list = "res")
  # run analysis on eeg and 3d eye data
  res <- mvjrqa(as.matrix(eye_3d),
                as.matrix(eeg),
                delay_1 = 1,
                embed_1 = 1,
                radius_1 = 10,
                delay_2 = 1,
                embed_2 = 1,
                radius_2 = 1,
                rescale = 0,
                normalize = 0,
                mindiagline = 2,
                minvertline = 2,
                tw = 1,
                whiteline = TRUE,
                recpt = FALSE,
                side = "both",
                metric = "euclidean",
                datatype = "continuous",
                setrec = TRUE,
                targetrec = rr_target)
  # store RR for individual RPs and JRP
  eeg_eye3d <- rbind(
    eeg_eye3d,
    data.frame(
      subject_id = data_files$subject_id[i],
      task_id = data_files$task_id[i],
      try = data_files$try[i],
      target_RR = rr_target,
      eye3d_RR1 = res[[1]]$RR,
      eye3d_RR2 = res[[2]]$RR,
      eye3d_JRR = res[[3]]$RR,
      eye3d_DET = res[[3]]$DET,
      eye3d_LAM = res[[3]]$LAM,
      eye3d_DENTR = res[[3]]$DENTR,
      eye3d_RAD1 = res[[1]]$radius,
      eye3d_RAD2 = res[[2]]$radius,
      eye3d_len1 = le1,
      eye3d_len2 = le2,
      eye3d_len3 = le3
    )
  )
  rm(list = "res")
  write.table(eeg_eye3d,
              file = paste0(eye3d_path, "/eye3d_", i, "_results.txt"),
              sep = ",", row.names = FALSE)
  }
  rm(list = c("eeg_eye3d", "eye_3d"))
  rm(list = c("eeg_eye2d", "eye_2d"))
  rm(list = "eeg")
}

## Parallelization
# Find out how many "workers" are available and prepare them for the task.
# Set num_workers to an appropriate value for your system at top of the script.
cl <-  makeCluster(num_workers)

# Create the same working environment for all of them.
clusterEvalQ(cl, {
  library("edf")
  source("R/mvjrqa.R")
  source("R/utils.R")
  NULL
})

clusterExport(cl, c("data_files", "eye2d_path", "eye3d_path"))

# Let them work
clusterApply(cl, i_vec, loop_data_fun)

# Send them home again (= shut down clusters).
stopCluster(cl)
