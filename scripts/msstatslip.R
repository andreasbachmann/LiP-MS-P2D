#!/usr/bin/env Rscript
# This script runs the MSstatsLiP pipeline on either test data (with '--test' flag)
# or real data (no flag)

# parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
test_mode <- "--test" %in% args

# set directory, input, output files and parameters
if (test_mode) {
  # use test data
  lip_file <- "data/test/test_lip_raw.csv.gz"
  trp_file <- "data/test/test_trp_raw.csv.gz"
  fasta_file <- "data/test/test.fasta"

  lip_summarized_output <- "test/results/test_summarized_lip_data.csv"
  model_output <- "test/results/test_adjusted_lip_model.csv"

  log_dir <- "test/logs"
  results_dir <- "test/results"

  cat("RUNNING IN TEST MODE\n")
  cat(rep("=", 30), "\n\n", sep = "")

} else {
  # use real data
  lip_file <- "data/combined_lip_raw.csv"
  trp_file <- "data/combined_trp_raw.csv"
  fasta_file <- "data/UP000005640_9606Human202411.fasta"

  lip_summarized_output <- "results/summarized_lip_data.csv"
  model_output <- "results/adjusted_lip_model.csv"

  log_dir <- "logs"
  results_dir <- "results"
}

# create directories and setup logs
if (!dir.exists(log_dir)) {
  dir.create(log_dir, recursive = TRUE)
}
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

log_file <- file.path(log_dir, paste0("msstatslip_run_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
sink(log_file, split = TRUE)

# start reading libraries
cat("Loading libraries...\n")
library(MSstatsLiP)
library(tidyverse)
library(data.table)
cat("...libraries loaded\n\n")

# start reading files
cat("Reading input files...\n")
cat("\tLiP file:", lip_file, "\n")
LiPRawData <- fread(lip_file)
cat("\tLiP data:", nrow(LiPRawData), "rows\n")

cat("\tTrP file:", trp_file, "\n")
TrPRawData <- fread(trp_file)
cat("\tTrP data:", nrow(TrPRawData), "rows\n\n")
cat("...input files read\n\n")

# run msstatslip pipeline
cat("MSstatsLiP Pipeline Starting\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(rep("=", 30), "\n\n", sep = "")

# run converter
MSstatsLiP_data <- SpectronauttoMSstatsLiPFormat(LiPRawData, fasta_file, TrPRawData)

# run summarization
MSstatsLiP_summarized <- dataSummarizationLiP(MSstatsLiP_data, MBimpute = FALSE, use_log_file = FALSE, append = FALSE)

# write replicate files
cat("\n\nWriting replicate-level data...\n")
fwrite(MSstatsLiP_summarized[["LiP"]]$ProteinLevelData, lip_summarized_output)
cat("\tSaved:", lip_summarized_output, "\n\n")

# run model
MSstatsLiP_model <- groupComparisonLiP(MSstatsLiP_summarized, fasta = fasta_file, use_log_file = FALSE, append = FALSE)

# write model
cat("Saving final results...\n")
fwrite(MSstatsLiP_model[["Adjusted.LiP.Model"]], model_output)
cat("\tSaved:", model_output, "\n\n")

cat("Pipeline Complete!\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(rep("=", 30), "\n")

sink()
cat("\nLog saved to:", log_file, "\n")
