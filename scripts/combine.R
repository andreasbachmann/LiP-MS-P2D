#!/usr/bin/env Rscript
# This script combines peptide-level p-values to domain-level p-values
# using Cauchy Combination Test and Empirical Brown's Method

# parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
test_mode <- "--test" %in% args
global_mode <- "--global" %in% args

# set mode suffix
mode_suffix <- ifelse(global_mode, "_global", "")

# set directory, input, output files and parameters
if (test_mode) {
  # use test data
  domain_annotations_file <- "data/test/test_domain_annotations.csv"
  msstatslip_model_file <- "test/results/test_adjusted_lip_model.csv"
  replicate_file <- "test/results/test_summarized_lip_data.csv"

  cauchy_output_file <- paste0("test/results/test_domain_level_results_cauchy", mode_suffix, ".csv")
  ebm_output_file <- paste0("test/results/test_domain_level_results_ebm", mode_suffix, ".csv")
  cauchy_dedup_file <- paste0("test/results/test_domain_level_results_cauchy_dedup", mode_suffix, ".csv")
  ebm_dedup_file <- paste0("test/results/test_domain_level_results_ebm_dedup", mode_suffix, ".csv")
  mapping_file <- paste0("test/results/test_peptide_domain_mapping", mode_suffix, ".csv")

  log_dir <- "test/logs"

  cat("\n")
  cat("RUNNING IN TEST MODE\n")
  cat(rep("=", 70), "\n", sep = "")

} else {
  # use real data
  msstatslip_model_file <- "results/adjusted_lip_model.csv"
  domain_annotations_file <- "data/domain_annotations_consolidated.csv"
  replicate_file <- "results/summarized_lip_data.csv"

  cauchy_output_file <- paste0("results/domain_level_results_cauchy", mode_suffix, ".csv")
  ebm_output_file <- paste0("results/domain_level_results_ebm", mode_suffix, ".csv")
  cauchy_dedup_file <- paste0("results/domain_level_results_cauchy_dedup", mode_suffix, ".csv")
  ebm_dedup_file <- paste0("results/domain_level_results_ebm_dedup", mode_suffix, ".csv")
  mapping_file <- paste0("results/peptide_domain_mapping.csv")

  log_dir <- "logs"
}

# setup logs
log_file <- file.path(log_dir, paste0("combineR_run", mode_suffix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
sink(log_file, split = TRUE)

# initiate start time
script_start <- Sys.time()

# print header
cat("\n")
if (global_mode) {
  cat("Domain-Level p-value Aggregation (GLOBAL MODE: by DomainName)\n")
} else {
  cat("Domain-Level p-value Aggregation (by DomainID)\n")
}
cat(rep("=", 70), "\n\n", sep = "")


# read libraries
cat("Loading required libraries...\n")
required_packages <- c("dplyr", "tidyr", "data.table", "EmpiricalBrownsMethod", "ACAT")

tryCatch(
  {
    suppressPackageStartupMessages({
      lapply(required_packages, function(pkg) {
        library(pkg, character.only = TRUE, quietly = TRUE)
      })
    })
    cat(" > libraries loaded successfully\n\n")
  },
  error = function(e) {
    cat("ERROR: Failed to load libraries\n")
    cat("Error message:", e$message, "\n")
    quit(status = 1)
  }
)

# read functions
source("R/aggregate_pvalues.R")
source("R/classify_trypticity.R")
source("R/deduplicate.R")
source("R/map_peptides_to_domains.R")
source("R/progress_bar.R")

# read files
cat("Reading files...\n")
peptide_results <- fread(msstatslip_model_file)
cat(sprintf("%-20s %25s\n", " > model as", msstatslip_model_file))
domain_annotations <- fread(domain_annotations_file)
cat(sprintf("%-20s %25s\n", " > domains as", domain_annotations_file))
replicate_data <- fread(replicate_file)
cat(sprintf("%-20s %25s\n\n", " > replicates as", replicate_file))

# run aggregation with appropriate mode
aggregate_by <- ifelse(global_mode, "domain_name", "domain_id")

result <- aggregate_pvalues(
  peptide_results,
  replicate_data,
  domain_annotations,
  aggregate_by = aggregate_by,
  save_mapping = TRUE,
  mapping_file = mapping_file
)

# deduplicate results
cat("\n")
cat("Deduplication\n")
cat(rep("=", 70), "\n\n", sep = "")

# deduplicate EBM
n_ebm_before <- nrow(result$ebm)
result$ebm_dedup <- deduplicate_results(result$ebm, global_mode = global_mode)
n_ebm_after <- nrow(result$ebm_dedup)
n_ebm_removed <- n_ebm_before - n_ebm_after
pct_ebm_removed <- 100 * n_ebm_removed / n_ebm_before

cat(sprintf("EBM deduplication:\n"))
cat(sprintf(" > before: %d\n", n_ebm_before))
cat(sprintf(" > after:  %d\n", n_ebm_after))
cat(sprintf(" > removed: %d (%.1f%%)\n\n", n_ebm_removed, pct_ebm_removed))

# deduplicate Cauchy
n_cauchy_before <- nrow(result$cauchy)
result$cauchy_dedup <- deduplicate_results(result$cauchy, global_mode = global_mode)
n_cauchy_after <- nrow(result$cauchy_dedup)
n_cauchy_removed <- n_cauchy_before - n_cauchy_after
pct_cauchy_removed <- 100 * n_cauchy_removed / n_cauchy_before

cat(sprintf("Cauchy deduplication:\n"))
cat(sprintf(" > before: %d\n", n_cauchy_before))
cat(sprintf(" > after:  %d\n", n_cauchy_after))
cat(sprintf(" > removed: %d (%.1f%%)\n\n\n", n_cauchy_removed, pct_cauchy_removed))

# save results
cat("Saving Results\n")
cat(rep("=", 70), "\n\n", sep = "")

# calculate runtime
total_elapsed <- as.numeric(difftime(Sys.time(), script_start, units = "secs"))
cat(sprintf("Total run time: %.1f seconds (%.1f minutes)\n\n", total_elapsed, total_elapsed / 60))

# save normal results
fwrite(result$cauchy, cauchy_output_file)
cat("Saved Cauchy results as", cauchy_output_file, "\n")
fwrite(result$ebm, ebm_output_file)
cat("Saved EBM results as", ebm_output_file, "\n")

# save deduplicated results
fwrite(result$cauchy_dedup, cauchy_dedup_file)
cat("Saved deduplicated Cauchy results as", cauchy_dedup_file, "\n")
fwrite(result$ebm_dedup, ebm_dedup_file)
cat("Saved deduplicated EBM results as", ebm_dedup_file, "\n")

# mapping file already saved inside aggregate_pvalues()
cat("Saved peptide-domain mapping as", mapping_file, "\n")

sink()
cat("Saved log as", log_file, "\n\n")