# process each direction with EBM
process_direction_ebm <- function(peptides, replicate_data, direction, 
											 verbose = FALSE, domain_id = NULL) {
	
	# which peptides are in the domain instance?
	requested_peptides <- peptides$FULL_PEPTIDE

	# which peptides that are in the domain instance are also in the replicate matrix?
	direction_replicates <- replicate_data %>%
		filter(FULL_PEPTIDE %in% requested_peptides)

	# only combine peptides that have replicates
	available_peptides <- direction_replicates$FULL_PEPTIDE
	sample_names <- setdiff(colnames(direction_replicates), "FULL_PEPTIDE")
	direction_matrix <- direction_replicates %>%
		select(-FULL_PEPTIDE) %>%
		as.matrix()

	# check for complete vs incomplete rows in matrix
	complete_rows <- complete.cases(direction_matrix)
	n_complete <- sum(complete_rows)
	n_incomplete <- sum(!complete_rows)

	incomplete_peptides <- NULL
	incomplete_peptide_nas <- NULL
	if (n_incomplete > 0) {
		incomplete_peptides <- available_peptides[!complete_rows]
		incomplete_matrix <- direction_matrix[!complete_rows, , drop = FALSE]
		na_counts_per_sample <- colSums(is.na(incomplete_matrix))
		incomplete_peptide_nas <- setNames(na_counts_per_sample, sample_names)
	}

	# only more than 2 peptides per direction can be combined
	if (n_complete < 2) {
		return(list(status = "less_than_two", n_incomplete = n_incomplete, incomplete_peptides = incomplete_peptides, incomplete_peptide_nas = incomplete_peptide_nas))
	}

	# filter to only complete peptides
	direction_matrix <- direction_matrix[complete_rows, , drop = FALSE]
	available_peptides <- available_peptides[complete_rows]
	peptides <- peptides %>% filter(FULL_PEPTIDE %in% available_peptides)

	# order pvalues to match matrix
	pvals_ordered <- peptides$adj.pvalue[match(available_peptides, peptides$FULL_PEPTIDE)] / 2

	# run EBM
	tryCatch({
		ebm_result <- empiricalBrownsMethod(
														data_matrix = direction_matrix,
														p_values = pvals_ordered,
														extra_info = TRUE)

		ebm_result$status <- "success"
		ebm_result$n_incomplete <- n_incomplete
		ebm_result$incomplete_peptides <- incomplete_peptides
		ebm_result$incomplete_peptide_nas <- incomplete_peptide_nas
		return(ebm_result)

	}, error = function(e) {
		return(list(status = "ebm_error", error_msg = e$message, n_incomplete = n_incomplete, incomplete_peptide_nas = incomplete_peptide_nas))
	})
}
