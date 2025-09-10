# FIXME: Remove the DRYRUN option before comitting

## ---------------------------
##
## Script name: 01_mvjrqa_analysis_empirical.R
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
## In order to run this script, you must first download the empirical data set
## by Shafei & Shadpour (2023) from Physionet.
##
## Link: https://doi.org/10.13026/kyjw-p786
##
## Note that this is a large data set with a size of 7.8 GB.
## 
## After downloading, place the contents in a directory called
## "empirical_data" that is a sub-directory of the directory where this
## R script is located.
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
## │   ├── RECORDS
##
## ---------------------------



# TODO: Make sure file passes lint test

# wipe workspace
rm(list = ls())

# load libraries and source functions
library(readr)
library(dplyr)
library(edf)
source("R/mvjrqa.R")
source("R/utils.R")
library(parallel)

# Individual files are written to results/eye2d and results/eye3d
eye2d_path <- "results/eye2d"
eye3d_path <- "results/eye3d"

# Utility function to create directories in the destination directory
# if they do note already exist
create_dir_if_not_present <- function(directory = NULL) {
  # Check that directory is a character variable
  if (!is.character(directory)) {
    stop("directory must be of type character")
  }
  # Check if a file with the same name exists, but is not a directory
  if (file.exists(directory) && !dir.exists(directory)) {
    stop(
      "A file: ",
      directory,
      " already exists, so a directory with the same name is not created."
    )
  }
  # Check if directory does not exist
  if (!dir.exists(directory)) {
    # Check if parent directory exists. If not, then create it
    dir_path <- unlist(strsplit(directory, "/"))
    path_length <- length(dir_path)
    # The parent directory is the path element up to the penultimate element
    # separated by "/".
    if (path_length > 1) {
      # Here basename() could be used instead
      parent_dir <- paste0(dir_path[1:path_length - 1], collapse="/")
      create_dir_if_not_present(parent_dir) # Recursively call the function
    }
    # Now the parent directory should exist, so `directory` can be created
    dir.create(directory)
  }
}

# Create directories folders "results/eye2d" and "results/eye3d"  if they do not exist
create_dir_if_not_present(eye2d_path)
create_dir_if_not_present(eye3d_path)

# Exclude the two EEG files and the two eye files that are not used.
eeg_exclude <- c("21_3_4.edf", "21_3_5.edf")
# Get the file names from the performance data file
perf <- read_csv("empirical_data/PerformanceScores.csv", quote = "'") |> 
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

# Create a vector of row index for parallellization
i_vec = 1:nrow(data_files) # vector defining the data files for loop_data_fun

### Function for parallelization
loop_data_fun = function(i, x) {
  DRYRUN <- TRUE # When set the mvjrqa will be skipped to allow testing the overall script function
  # load and sort eeg data
  # TODO: Add directory to path for eeg_data[] file names. Should be done with the paste0()
  eeg <- read.edf(paste0("empirical_data/EEG_FLS/", data_files$eeg_data[i]), read.annotations = TRUE, header.only = FALSE) # load eeg data
  eeg <- data.frame(eeg$signal) # convert to data frame
  eeg <- eeg[seq(1,dim(eeg)[1], by = 10),seq(1,dim(eeg)[2], by = 2)] # remove t-variables and downsample to 50 Hz
  drops <- c("EEGHEOGRCPz.data", "EEGHEOGLCPz.data", "EEGVEOGUCPz.data","EEGVEOGLCPz.data") # list EOGs to drop from data frame
  eeg <- eeg[ ,!(names(eeg) %in% drops)] # drop EOGs from data frame
  rm(list = "drops")
  
  # load eye data
  df_eye <- read.csv(paste0("empirical_data/EYE_FLS/",data_files$eye_data[i])) # load eye movement data
  eye_2d <- data.frame(df_eye[,1:2]) # extract 2d eye gaze
  eye_3d <- data.frame(df_eye[,3:5]) # extract 3d gaze
  
  # cut off superfluous data points
  if(dim(eeg)[1] > dim(df_eye)[1]) {
    eeg <- eeg[1:dim(eye_2d)[1],] # cut eeg record
  } else if (dim(eeg)[1] < dim(df_eye)[1]) {
    eye_2d <- eye_2d[1:dim(eeg)[1],]# cut eye records
    eye_3d <- eye_3d[1:dim(eeg)[1],]# cut eye records
  } else {
    # do nothing
  }
  le1 <- dim(eeg)[1]
  rm(list = "df_eye")
  
  # remove NAs from all data sets
  naCount <- rowSums(is.na(eeg)) + rowSums(is.na(eye_2d)) + rowSums(is.na(eye_3d))
  eeg <- eeg[naCount==0,]
  eye_2d <- eye_2d[naCount==0,]
  eye_3d <- eye_3d[naCount==0,]
  rm(list = "naCount")
  le2 <- dim(eeg)[1]
  
  # Calculate difference values for each data set
  eeg <- data.frame(lapply(eeg, diff))
  eye_2d <- data.frame(lapply(eye_2d, diff))
  eye_3d <- data.frame(lapply(eye_3d, diff))
  
  # Calculate means and SDs for EEG data
  stats <- sapply(eeg, function(x) c(mean = mean(x), sd = sd(x)))
  
  # remove variables with 0 SD
  eeg <- eeg[,stats[2,] > 0]
  stats <- stats[,stats[2,] > 0]
  
  # remove data outside +/-3 SD
  outside_bounds <- as.data.frame(lapply(names(eeg), function(var) {
    # Accessing eeg directly and stats by row and column names
    abs(eeg[[var]] - stats["mean", var]) >= 3 * stats["sd", var]
  }))
  
  # create variable that codes for rows with no data point outside of 3SD
  SD3count <- rowSums(outside_bounds)
  
  # remove data accordingly
  eeg <- eeg[SD3count==0,]
  eye_2d <- eye_2d[SD3count==0,]
  eye_3d <- eye_3d[SD3count==0,]
  rm(list = c("stats","outside_bounds", "SD3count"))
  
  # length of final data set
  le3 <- dim(eeg)[1]
  
  # run analysis on eeg and 2d eye data
  if (DRYRUN) {
    # If DRYRUN is TRUE use test data
    eeg_eye2d <- data.frame(
      subject_id = data_files$subject_id[i],
      task_id = data_files$task_id[i],
      try = data_files$try[i],
      eye2d_RR1 = 1,
      eye2d_RR2 = 2,
      eye2d_JRR = 3,
      eye2d_DET = 4,
      eye2d_LAM = 5,
      eye2d_DENTR = 6,
      eye2d_RAD1 = 7,
      eye2d_RAD2 = 8,
      eye2d_len1 = 9,
      eye2d_len2 = 10,
      eye2d_len3 = 11
    )
  } else {
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
                  targetrec = 5)
    
    # store RR for individual RPs and JRP
    # TODO: Give the data frame column names
    # TODO: Add subject_id, task_id, try
    eeg_eye2d <- data.frame(
      subject_id = data_files$subject_id[i],
      task_id = data_files$task_id[i],
      try = data_files$try[i],
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
  }
  write.table(eeg_eye2d, file = paste0(eye2d_path, "/eye2d_", i, "_results.txt"), sep = ",", row.names = FALSE)
  rm(list = c("res", "eeg_eye2d", "eye_2d"))
  
  # run analysis on eeg and 3d eye data
  if (DRYRUN) {
    # If DRYRUN is TRUE use test data
    eeg_eye3d <- data.frame(
      subject_id = data_files$subject_id[i],
      task_id = data_files$task_id[i],
      try = data_files$try[i],
      eye3d_RR1 = 1,
      eye3d_RR2 = 2,
      eye3d_JRR = 3,
      eye3d_DET = 4,
      eye3d_LAM = 5,
      eye3d_DENTR = 6,
      eye3d_RAD1 = 7,
      eye3d_RAD2 = 8,
      eye3d_len1 = 9,
      eye3d_len2 = 10,
      eye3d_len3 = 11
    )
  } else {
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
                  targetrec = 5)
    
    # store RR for individual RPs and JRP
    # TODO: Give the data frame column names
    # TODO: Add subject_id, task_id, try
    eeg_eye3d <- data.frame(
      subject_id = data_files$subject_id[i],
      task_id = data_files$task_id[i],
      try = data_files$try[i],
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
  }
  write.table(eeg_eye3d, file = paste0(eye3d_path, "/eye3d_", i, "_results.txt"), sep = ",", row.names = FALSE)
  rm(list = c("res", "eeg_eye3d", "eye_3d"))
  rm(list = "eeg")
}

## Parallelization
# Find out how many "workers" are available and prepare them for the task.
cl = makeCluster(6)

# Create the same working environment for all of them.
clusterEvalQ(cl, {
  library("edf")
  source("R/mvjrqa.R")
  source("R/utils.R")})

clusterExport(cl, c("data_files", "eye2d_path", "eye3d_path"))

# Let them work
clusterApply(cl, i_vec, loop_data_fun)

# Send them home again (= shut down clusters).
stopCluster(cl)