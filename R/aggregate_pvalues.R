# main function aggregating peptide-level p-values to domain-level with Cauchy and EBM
aggregate_pvalues <- function(peptide_results, replicate_data, domain_annotations,
                              overlap_cutoff = 0.5, pval_cutoff = 0.05,
                              verbose = FALSE, save_mapping = FALSE, mapping_file = "",
                              aggregate_by = "domain_id") {
  # validate aggregate_by parameter
  if (!aggregate_by %in% c("domain_id", "domain_name")) {
    stop("aggregate_by must be 'domain_id' or 'domain_name'")
  }

  # set grouping column based on mode
  group_col <- ifelse(aggregate_by == "domain_id", "DomainID", "DomainName")

  ######################
  # TRYPTICITY FILTERING
  cat("Classifying trypticity...\n")
  peptide_results <- classify_trypticity(peptide_results)

  n_total <- nrow(peptide_results)
  n_semi <- sum(peptide_results$trypticity == "semi_tryptic")
  n_fully <- sum(peptide_results$trypticity == "fully_tryptic")

  cat(sprintf(" > Total peptides: %d\n", n_total))
  cat(sprintf(" > Semi-tryptic: %d (%.1f%%)\n", n_semi, 100 * n_semi / n_total))
  cat(sprintf(" > Fully-tryptic: %d (%.1f%%)\n", n_fully, 100 * n_fully / n_total))

  peptide_results <- peptide_results %>%
    filter(trypticity == "semi_tryptic")

  cat(sprintf(" > Kept %d semi-tryptic peptides for analysis\n\n", nrow(peptide_results)))

  ################
  # DOMAIN MAPPING
  cat("Mapping peptides to domains...\n")
  cat(sprintf(" > Overlap cutoff: %d%%\n", overlap_cutoff * 100))

  peptide_domain_mapping <- map_peptides_to_domains(peptide_results, domain_annotations, overlap_cutoff)

  mapped_peptides <- unique(peptide_domain_mapping$FULL_PEPTIDE)
  n_mapped <- length(mapped_peptides)
  n_total_pep <- n_distinct(peptide_results$FULL_PEPTIDE)
  n_domains_mapped <- n_distinct(peptide_domain_mapping$DomainID)

  pep_per_domain <- peptide_domain_mapping %>%
    group_by(across(all_of(group_col))) %>%
    summarise(n = n(), .groups = "drop")

  cat(sprintf(
    " > Peptides per %s: mean=%.1f, median=%d, max=%d\n",
    group_col, mean(pep_per_domain$n), median(pep_per_domain$n), max(pep_per_domain$n)
  ))
  cat(sprintf(
    " > Mapped %d / %d peptides (%.1f%% unmapped)\n",
    n_mapped, n_total_pep, 100 * (n_total_pep - n_mapped) / n_total_pep
  ))
  cat(sprintf(" > To %d domain instances\n", n_domains_mapped))
  cat(sprintf(" > Total peptide-domain links: %d\n", nrow(peptide_domain_mapping)))

  # optionally save mapping
  if (save_mapping && mapping_file != "") {
    peptide_details <- peptide_domain_mapping %>%
      inner_join(
        peptide_results %>%
          select(FULL_PEPTIDE, Label, log2FC, adj.pvalue, Tvalue, DF, trypticity),
        by = "FULL_PEPTIDE",
        relationship = "many-to-many"
      ) %>%
      filter(
        !is.na(Tvalue), is.finite(Tvalue),
        !is.na(DF), is.finite(DF),
        DF > 0
      ) %>%
      mutate(
        p_right = pt(Tvalue, df = DF, lower.tail = FALSE),
        p_left = pt(Tvalue, df = DF, lower.tail = TRUE)
      ) %>%
      arrange(across(all_of(group_col)), Label)

    fwrite(peptide_details, file = mapping_file)
    cat(sprintf(" > Saved mapping to %s\n\n", mapping_file))
  }

  ################
  # CONTRAST SETUP
  contrasts <- unique(peptide_results$Label)
  cat(sprintf("Found %d contrasts...\n\n", length(contrasts)))

  results_ebm <- list()
  results_cauchy <- list()

  ###############
  # CONTRAST LOOP
  for (contrast_idx in seq_along(contrasts)) {
    contrast <- contrasts[contrast_idx]

    cat(sprintf("Contrast %d/%d: %s\n", contrast_idx, length(contrasts), contrast))
    cat(rep("=", 70), "\n", sep = "")

    conditions <- trimws(strsplit(contrast, " vs ")[[1]])

    #####################
    # FILTER & ONE-TAILED
    peptide_pvals <- peptide_results %>%
      filter(Label == contrast) %>%
      filter(FULL_PEPTIDE %in% mapped_peptides) %>%
      filter(
        !is.na(log2FC), is.finite(log2FC),
        !is.na(Tvalue), is.finite(Tvalue),
        !is.na(DF), is.finite(DF),
        DF > 0
      ) %>%
      mutate(
        p_right = pt(Tvalue, df = DF, lower.tail = FALSE),
        p_left = pt(Tvalue, df = DF, lower.tail = TRUE)
      ) %>%
      select(FULL_PEPTIDE, log2FC, adj.pvalue, p_right, p_left, Tvalue, DF, trypticity)

    cat(sprintf(" > Valid peptides: %d\n", nrow(peptide_pvals)))

    # join with domain mapping
    contrast_peptides <- peptide_domain_mapping %>%
      inner_join(peptide_pvals, by = "FULL_PEPTIDE")

    group_ids <- unique(contrast_peptides[[group_col]])
    n_groups <- length(group_ids)
    cat(sprintf(" > Domains to process: %d\n", n_groups))

    ##################
    # REPLICATE MATRIX
    contrast_replicates <- replicate_data %>%
      filter(GROUP %in% conditions) %>%
      filter(FULL_PEPTIDE %in% peptide_pvals$FULL_PEPTIDE) %>%
      mutate(SAMPLE_ID = paste(GROUP, SUBJECT, sep = "_")) %>%
      select(FULL_PEPTIDE, SAMPLE_ID, LogIntensities) %>%
      pivot_wider(names_from = SAMPLE_ID, values_from = LogIntensities, id_cols = FULL_PEPTIDE)

    cat(sprintf(
      " > Replicate matrix: %d peptides x %d samples\n\n",
      nrow(contrast_replicates), ncol(contrast_replicates) - 1
    ))

    # tracking
    stats <- list(
      n_single = 0,
      n_multi = 0,
      ebm_attempted = 0,
      ebm_success = 0,
      ebm_failed = 0,
      ebm_full = 0,
      ebm_reduced = 0,
      cauchy_attempted = 0,
      cauchy_success = 0,
      cauchy_failed = 0
    )

    contrast_results_ebm <- list()
    contrast_results_cauchy <- list()

    cat("Processing domains...\n")

    # initialize progress bar
    cat(sprintf(" > %5d/%5d - %s", 0, n_groups, create_progress_bar(0)))

    #############
    # DOMAIN LOOP
    for (i in seq_along(group_ids)) {
      group_id <- group_ids[i]

      # print progress
      if (runif(1) < 0.1 || i == n_groups || i == 1) {
        cat(sprintf("\r > %5d/%5d - %s", i, n_groups, create_progress_bar(100 * i / n_groups)))
      }

      # get peptides for this domain
      group_peptides <- contrast_peptides %>%
        filter(.data[[group_col]] == group_id)

      n_pep <- nrow(group_peptides)
      if (n_pep == 0) next

      # metadata
      if (aggregate_by == "domain_name") {
        domain_name <- group_id
        protein_name <- paste(unique(group_peptides$ProteinName), collapse = ";")
        domain_id <- NA
        n_proteins <- n_distinct(group_peptides$ProteinName)
        n_instances <- n_distinct(group_peptides$DomainID)
      } else {
        domain_name <- group_peptides$DomainName[1]
        protein_name <- group_peptides$ProteinName[1]
        domain_id <- group_id
        n_proteins <- NA
        n_instances <- NA
      }

      # summary stats
      mean_fc <- mean(group_peptides$log2FC)
      median_fc <- median(group_peptides$log2FC)

      # base result template
      make_result <- function(p_right, p_left, method, n_used = n_pep) {
        res <- data.frame(
          DomainID = domain_id,
          DomainName = domain_name,
          ProteinName = protein_name,
          Label = contrast,
          n_peptides = n_pep,
          n_peptides_used = n_used,
          mean_log2FC = mean_fc,
          median_log2FC = median_fc,
          p_right = p_right,
          p_left = p_left,
          method = method,
          stringsAsFactors = FALSE
        )
        if (aggregate_by == "domain_name") {
          res$n_proteins <- n_proteins
          res$n_domain_instances <- n_instances
        }
        return(res)
      }

      ################
      # SINGLE PEPTIDE
      if (n_pep == 1) {
        stats$n_single <- stats$n_single + 1

        single_result <- make_result(
          p_right = group_peptides$p_right[1],
          p_left = group_peptides$p_left[1],
          method = "single_peptide"
        )

        contrast_results_ebm[[group_id]] <- single_result
        contrast_results_cauchy[[group_id]] <- single_result
        next
      }

      ###############
      # MULTI PEPTIDE
      stats$n_multi <- stats$n_multi + 1

      ########
      # CAUCHY
      stats$cauchy_attempted <- stats$cauchy_attempted + 1

      tryCatch(
        {
          cauchy_right <- ACAT(group_peptides$p_right)
          cauchy_left <- ACAT(group_peptides$p_left)

          contrast_results_cauchy[[group_id]] <- make_result(
            p_right = cauchy_right,
            p_left = cauchy_left,
            method = "cauchy"
          )
          stats$cauchy_success <- stats$cauchy_success + 1
        },
        error = function(e) {
          stats$cauchy_failed <<- stats$cauchy_failed + 1
          if (verbose) cat(sprintf("\n > Cauchy failed for %s: %s\n", group_id, e$message))
        }
      )

      #####
      # EBM
      stats$ebm_attempted <- stats$ebm_attempted + 1

      # get replicate data for this domain's peptides
      domain_replicates <- contrast_replicates %>%
        filter(FULL_PEPTIDE %in% group_peptides$FULL_PEPTIDE)

      if (nrow(domain_replicates) < 2) {
        stats$ebm_failed <- stats$ebm_failed + 1
        next
      }

      # build matrix
      rep_matrix <- domain_replicates %>%
        select(-FULL_PEPTIDE) %>%
        as.matrix()

      available_peptides <- domain_replicates$FULL_PEPTIDE

      # filter to complete cases
      complete_rows <- complete.cases(rep_matrix)
      if (sum(complete_rows) < 2) {
        stats$ebm_failed <- stats$ebm_failed + 1
        next
      }

      rep_matrix <- rep_matrix[complete_rows, , drop = FALSE]
      available_peptides <- available_peptides[complete_rows]
      n_used_ebm <- length(available_peptides)

      # get p-values in correct order
      pep_subset <- group_peptides %>% filter(FULL_PEPTIDE %in% available_peptides)
      idx <- match(available_peptides, pep_subset$FULL_PEPTIDE)
      p_right_ordered <- pep_subset$p_right[idx]
      p_left_ordered <- pep_subset$p_left[idx]


      # EBM for both tails
      ebm_right <- NULL
      ebm_left <- NULL

      tryCatch(
        {
          ebm_result_right <- empiricalBrownsMethod(
            data_matrix = rep_matrix,
            p_values = p_right_ordered,
            extra_info = TRUE
          )
          ebm_right <- ebm_result_right$P_test
        },
        error = function(e) {
          if (verbose) cat(sprintf("\n > EBM right failed for %s: %s\n", group_id, e$message))
        }
      )

      tryCatch(
        {
          ebm_result_left <- empiricalBrownsMethod(
            data_matrix = rep_matrix,
            p_values = p_left_ordered,
            extra_info = TRUE
          )
          ebm_left <- ebm_result_left$P_test
        },
        error = function(e) {
          if (verbose) cat(sprintf("\n > EBM left failed for %s: %s\n", group_id, e$message))
        }
      )

      # only store if both succeeded
      if (!is.null(ebm_right) && !is.null(ebm_left)) {
        contrast_results_ebm[[group_id]] <- make_result(
          p_right = ebm_right,
          p_left = ebm_left,
          method = "EBM",
          n_used = n_used_ebm
        )
        stats$ebm_success <- stats$ebm_success + 1
        if (n_used_ebm < n_pep) {
          stats$ebm_reduced <- stats$ebm_reduced + 1
        } else {
          stats$ebm_full <- stats$ebm_full + 1
        }
      } else {
        stats$ebm_failed <- stats$ebm_failed + 1
      }
    }

    cat("\n\n")

    ##################
    # CONTRAST SUMMARY
    cat("Domain breakdown:\n")
    cat(sprintf(" > Single peptide: %d\n", stats$n_single))
    cat(sprintf(" > Multi peptide: %d\n", stats$n_multi))

    cat(sprintf(
      "\nEBM: %d/%d successful (%.1f%%)\n",
      stats$ebm_success, stats$ebm_attempted,
      100 * stats$ebm_success / max(stats$ebm_attempted, 1)
    ))

    cat(sprintf(" > domains failed (< 2 complete):    %d (%.1f%%)\n", stats$ebm_failed, 100 * stats$ebm_failed / max(stats$ebm_attempted, 1)))

    cat(sprintf(
      " > domains using all peptides:  %d (%.1f%%)\n",
      stats$ebm_full, 100 * stats$ebm_full / max(stats$ebm_success, 1)
    ))
    cat(sprintf(
      " > domains dropping some peptides:   %d (%.1f%%)\n",
      stats$ebm_reduced, 100 * stats$ebm_reduced / max(stats$ebm_success, 1)
    ))


    cat(sprintf(
      "Cauchy: %d/%d successful (%.1f%%)\n\n",
      stats$cauchy_success, stats$cauchy_attempted,
      100 * stats$cauchy_success / max(stats$cauchy_attempted, 1)
    ))


    # store results
    if (length(contrast_results_ebm) > 0) {
      results_ebm[[contrast]] <- bind_rows(contrast_results_ebm)
    }
    if (length(contrast_results_cauchy) > 0) {
      results_cauchy[[contrast]] <- bind_rows(contrast_results_cauchy)
    }
  }

  ###############
  # COMBINE & FDR
  cat("Combining results and computing FDR...\n")
  cat(rep("=", 70), "\n", sep = "")

  final_ebm <- bind_rows(results_ebm)
  final_cauchy <- bind_rows(results_cauchy)

  # BH correction for each tail separately within each contrast
  if (nrow(final_ebm) > 0) {
    final_ebm <- final_ebm %>%
      group_by(Label) %>%
      mutate(
        adj_p_right = p.adjust(p_right, method = "BH"),
        adj_p_left = p.adjust(p_left, method = "BH")
      ) %>%
      ungroup() %>%
      arrange(pmin(p_right, p_left))

    n_sig_right <- sum(final_ebm$adj_p_right < pval_cutoff, na.rm = TRUE)
    n_sig_left <- sum(final_ebm$adj_p_left < pval_cutoff, na.rm = TRUE)
    cat(sprintf(
      " > EBM: %d domains, %d significant right (FDR < %.2f), %d significant left\n",
      nrow(final_ebm), n_sig_right, pval_cutoff, n_sig_left
    ))
  }

  if (nrow(final_cauchy) > 0) {
    final_cauchy <- final_cauchy %>%
      group_by(Label) %>%
      mutate(
        adj_p_right = p.adjust(p_right, method = "BH"),
        adj_p_left = p.adjust(p_left, method = "BH")
      ) %>%
      ungroup() %>%
      arrange(pmin(p_right, p_left))

    n_sig_right <- sum(final_cauchy$adj_p_right < pval_cutoff, na.rm = TRUE)
    n_sig_left <- sum(final_cauchy$adj_p_left < pval_cutoff, na.rm = TRUE)
    cat(sprintf(
      " > Cauchy: %d domains, %d significant right (FDR < %.2f), %d significant left\n",
      nrow(final_cauchy), n_sig_right, pval_cutoff, n_sig_left
    ))
  }

  cat("\nPeptide to Domain Aggregation completed!\n")

  return(list(
    ebm = final_ebm,
    cauchy = final_cauchy
  ))
}