# main function aggregating peptide-level p-values to domain-level with Cauchy and EBM
aggregate_pvalues <- function(peptide_results, replicate_data, domain_annotations,
                              overlap_cutoff = 0.5, pval_cutoff = 0.05, fc_cutoff = 1,
                              verbose = FALSE, impute = FALSE, save_mapping = FALSE, mapping_file = "") {
  
  # trypticity classification
  cat("Classifying trypticity...\n")
  peptide_results <- classify_trypticity(peptide_results)

  n_total <- nrow(peptide_results)
  n_semi <- sum(peptide_results$trypticity == "semi_tryptic")
  n_fully <- sum(peptide_results$trypticity == "fully_tryptic")

  # output trypticity stats
  cat(sprintf(" > Total peptides: %d\n", n_total))
  cat(sprintf(" > Semi-tryptic: %d (%.1f%%)\n", n_semi, 100 * n_semi / n_total))
  cat(sprintf(" > Fully-tryptic: %d (%.1f%%)\n", n_fully, 100 * n_fully / n_total))

  # filter to semi-tryptic only
  peptide_results <- peptide_results %>%
    filter(trypticity == "semi_tryptic")

  cat(sprintf(" > Kept %d semi-tryptic peptides for analysis\n\n", nrow(peptide_results)))


  # domain mapping
  cat("Mapping peptides to domains...\n")
  cat(sprintf(" > overlap cutoff: %d%%\n", overlap_cutoff * 100))

  # run mapping
  peptide_domain_mapping <- map_peptides_to_domains(peptide_results, domain_annotations, overlap_cutoff)

  # key metrics
  mapped_peptides <- unique(peptide_domain_mapping$FULL_PEPTIDE)
  n_mapped <- length(mapped_peptides)
  n_total <- n_distinct(peptide_results$FULL_PEPTIDE)
  n_domains_mapped <- n_distinct(peptide_domain_mapping$DomainID)

  # peptides per domain distribution
  pep_per_domain <- peptide_domain_mapping %>%
    count(DomainID)

  # output mapping stats
  cat(sprintf(
    " > peptides per domain: mean=%.1f, median=%d, max=%d\n",
    mean(pep_per_domain$n), median(pep_per_domain$n), max(pep_per_domain$n)
  ))
  cat(sprintf(
    " > mapped %d / %d peptides (%.1f%% unmapped)\n",
    n_mapped, n_total, 100 * (n_total - n_mapped) / n_total
  ))
  cat(sprintf(
    " > to %d / %d domains (%.1f%% without coverage)\n",
    n_domains_mapped, nrow(domain_annotations),
    100 * (nrow(domain_annotations) - n_domains_mapped) / nrow(domain_annotations)
  ))
  cat(sprintf(" > total peptide-domain links: %d\n", nrow(peptide_domain_mapping)))

  if (save_mapping) {
    peptide_details <- peptide_domain_mapping %>%
      inner_join(
        peptide_results %>%
          select(FULL_PEPTIDE, Label, log2FC, adj.pvalue, trypticity),
        by = "FULL_PEPTIDE",
        relationship = "many-to-many"
      ) %>%
      mutate(Direction = ifelse(log2FC > 0, "UP", "DOWN")) %>%
      select(
        DomainID, DomainName, ProteinName, Label, Direction,
        FULL_PEPTIDE, PeptideSequence,
        log2FC, adj.pvalue,
        trypticity, StartPos, EndPos
      ) %>%
      arrange(DomainID, Label, Direction)
  
    fwrite(peptide_details, file = mapping_file)
    cat(sprintf(" > saved %d peptide x domain x contrast records\n\n", nrow(peptide_details)))
  }

  # get contrast info
  contrasts <- unique(peptide_results$Label)
  cat(sprintf("Found %d contrasts...\n\n", length(contrasts)))

  # global tracking & results lists
  global_cauchy_stats <- list()
  global_ebm_stats <- list()
  global_domain_stats <- list()
  results_by_contrast_cauchy <- list()
  results_by_contrast_ebm <- list()

  #################
  # CONTRAST LOOP #
  #################
  for (contrast_idx in seq_along(contrasts)) {
    contrast <- contrasts[contrast_idx]

    # contrast header
    cat(sprintf("Contrast %d/%d: %s\n", contrast_idx, length(contrasts), contrast))
    cat(rep("=", 70), "\n\n", sep = "")

    # parse current contrast
    conditions <- trimws(strsplit(contrast, " vs ")[[1]])

    # get p-values from all available peptides
    cat("Filtering Data...\n")
    peptide_pvals <- peptide_results %>%
      filter(Label == contrast) %>%
      filter(FULL_PEPTIDE %in% mapped_peptides) %>%
      filter(!is.na(log2FC), is.finite(log2FC), !is.na(adj.pvalue)) %>%
      select(FULL_PEPTIDE, log2FC, adj.pvalue, trypticity)

    valid_peptides <- unique(peptide_pvals$FULL_PEPTIDE)
    invalid_peptides <- peptide_results %>%
      filter(Label == contrast) %>%
      filter(FULL_PEPTIDE %in% mapped_peptides) %>%
      filter(is.na(log2FC) | !is.finite(log2FC) | is.na(adj.pvalue))

    # metrics; invalid data (NA, infinite)
    pct_invalid <- 100 * nrow(invalid_peptides) / n_mapped
    pct_valid <- 100 * nrow(peptide_pvals) / n_mapped
    cat(sprintf("%-42s %27s\n", " > excluded invalid peptides:", sprintf("(%.1f%%) %d / %d", pct_invalid, nrow(invalid_peptides), n_mapped)))
    cat(sprintf("%-42s %27s\n", " > valid peptides:", sprintf("(%.1f%%) %d / %d", pct_valid, nrow(peptide_pvals), n_mapped)))

    # join pvalues with domain mapping
    contrast_peptides <- peptide_domain_mapping %>%
      inner_join(peptide_pvals, by = "FULL_PEPTIDE")
    
    # metrics; get unique domains in current contrast, total peptide x domain matchings, (any domains lost, cause no pvalue?)
    domain_ids <- unique(contrast_peptides$DomainID)
    n_domains <- length(domain_ids)
    domains_lost_no_pvals <- n_domains_mapped - n_domains
    if (domains_lost_no_pvals > 0) {
      cat(sprintf("%-42s %27s\n", " > domains lost:", domains_lost_no_pvals))
    }
    cat(sprintf("%-42s %27s\n", " > domains to process:", n_domains))

    # prepare replicate data
    contrast_replicates <- replicate_data %>%
      filter(GROUP %in% c(conditions[1], conditions[2])) %>%
      filter(FULL_PEPTIDE %in% valid_peptides) %>%
      mutate(
        SAMPLE_ID = paste(GROUP, SUBJECT, sep = "_"),
        LogIntensities = as.numeric(LogIntensities)
      ) %>%
      select(FULL_PEPTIDE, SAMPLE_ID, LogIntensities) %>%
      pivot_wider(
        names_from = SAMPLE_ID,
        values_from = LogIntensities,
        id_cols = FULL_PEPTIDE
      )

    # metrics; replicate matrix dimensions & completion
    n_samples <- ncol(contrast_replicates) - 1
    n_peptides_matrix <- nrow(contrast_replicates)
    cat(sprintf(
      "%-42s %27s\n", " > replicate matrix dimensions:",
      sprintf("%d peptides x %d samples", n_peptides_matrix, n_samples)
    ))
    n_complete_peptides <- sum(complete.cases(contrast_replicates[, -1]))
    pct_complete <- 100 * n_complete_peptides / n_peptides_matrix
    cat(sprintf(
      "%-42s %27s\n\n", " > complete peptides in matrix:",
      sprintf("(%.1f%%) %d / %d", pct_complete, n_complete_peptides, n_peptides_matrix)
    ))

    # domain processing preparation
    cat("Processing Domains...\n")
    results_by_domain_cauchy <- list()
    results_by_domain_ebm <- list()

    domain_stats <- list(
      n_single_peptide = 0,
      n_multi_peptide = 0,
      n_up_only = 0,
      n_down_only = 0,
      n_mixed = 0
    )

    direction_stats <- list(
      n_total_instances = 0,
      n_single_peptide = 0,
      n_multi_peptide = 0
    )

    # tracking for cauchy method
    cauchy_stats <- list(
      n_attempted = 0,
      n_success = 0,
      n_failed = 0,
      n_single = 0
    )

    # tracking for ebm method
    ebm_stats <- list(
      n_attempted = 0,
      n_success = 0,
      n_failed_less_than_two = 0,
      n_failed_ebm_error = 0,
      n_peptides_removed = 0,
      unique_incomplete_peptides = character(),
      sample_na_counts = list()
    )

    # initialize progress bar
    cat(sprintf(" > %5d/%5d - %s", 0, n_domains, create_progress_bar(0)))

    ###############
    # DOMAIN LOOP #
    ###############

    for (i in seq_along(domain_ids)) {
      domain_id <- domain_ids[i]

      # print progress
      if (runif(1) < 0.1 || i == n_domains || i == 1) {
        cat(sprintf("\r > %5d/%5d - %s", i, n_domains, create_progress_bar(100 * i / n_domains)))
      }

      # get domain peptides
      domain_peptides <- contrast_peptides %>%
        filter(DomainID == domain_id)

      # peptide count in domain
      n_peptides <- nrow(domain_peptides)

      # how many peptides?
      # no valid peptides in domain
      if (n_peptides == 0) {
        cat("No peptide in domain... This shouldn't happen.")
        next
      }

      # single peptide in domain
      if (n_peptides == 1) {
        domain_stats$n_single_peptide <- domain_stats$n_single_peptide + 1
        direction_stats$n_total_instances <- direction_stats$n_total_instances + 1
        direction_stats$n_single_peptide <- direction_stats$n_single_peptide + 1

        one_tailed_pval <- domain_peptides$adj.pvalue[1] / 2
        direction <- ifelse(domain_peptides$log2FC[1] > 0, "UP", "DOWN")

        # same result for both methods when single peptide
        single_result <- data.frame(
          DomainID = domain_id,
          DomainName = domain_peptides$DomainName[1],
          ProteinName = domain_peptides$ProteinName[1],
          Label = contrast,
          Direction = direction,
          n_peptides = 1,
          combined_pvalue = one_tailed_pval,
          mean_log2FC = domain_peptides$log2FC[1],
          median_log2FC = domain_peptides$log2FC[1],
          method = "single_peptide",
          stringsAsFactors = FALSE
        )

        results_by_domain_ebm[[paste0(domain_id, "_", direction)]] <- single_result
        results_by_domain_cauchy[[paste0(domain_id, "_", direction)]] <- single_result
        next
      }

      # multi-peptide
      domain_stats$n_multi_peptide <- domain_stats$n_multi_peptide + 1

      # split peptides by direction
      up_peptides <- domain_peptides %>% filter(log2FC > 0)
      down_peptides <- domain_peptides %>% filter(log2FC < 0)
      n_up <- nrow(up_peptides)
      n_down <- nrow(down_peptides)

      # track direction type
      if (n_up > 0 && n_down == 0) {
        domain_stats$n_up_only <- domain_stats$n_up_only + 1
      } else if (n_down > 0 && n_up == 0) {
        domain_stats$n_down_only <- domain_stats$n_down_only + 1
      } else if (n_up > 0 && n_down > 0) {
        domain_stats$n_mixed <- domain_stats$n_mixed + 1
      }

      # process UP direction #
      ########################
      if (n_up > 0) {
        direction_stats$n_total_instances <- direction_stats$n_total_instances + 1

        # single peptide UP
        if (n_up == 1) {
          direction_stats$n_single_peptide <- direction_stats$n_single_peptide + 1

          one_tailed_pval <- up_peptides$adj.pvalue[1] / 2
          single_up_result <- data.frame(
            DomainID = domain_id,
            DomainName = up_peptides$DomainName[1],
            ProteinName = up_peptides$ProteinName[1],
            Label = contrast,
            Direction = "UP",
            n_peptides = 1,
            combined_pvalue = one_tailed_pval,
            mean_log2FC = up_peptides$log2FC[1],
            median_log2FC = up_peptides$log2FC[1],
            method = "single_peptide",
            stringsAsFactors = FALSE
          )
          results_by_domain_ebm[[paste0(domain_id, "_UP")]] <- single_up_result
          results_by_domain_cauchy[[paste0(domain_id, "_UP")]] <- single_up_result

          # multi-peptide UP
        } else if (n_up >= 2) {
          direction_stats$n_multi_peptide <- direction_stats$n_multi_peptide + 1

          # EBM method
          ebm_stats$n_attempted <- ebm_stats$n_attempted + 1
          ebm_result_up <- process_direction_ebm(
            peptides = up_peptides,
            replicate_data = contrast_replicates,
            direction = "UP",
            verbose = verbose,
            domain_id = domain_id
          )

          if (!is.null(ebm_result_up$status)) {
            if (ebm_result_up$n_incomplete > 0) {
              ebm_stats$n_peptides_removed <- ebm_stats$n_peptides_removed + ebm_result_up$n_incomplete

              # track unique incomplete peptides
              if (!is.null(ebm_result_up$incomplete_peptides)) {
                ebm_stats$unique_incomplete_peptides <- c(
                  ebm_stats$unique_incomplete_peptides,
                  ebm_result_up$incomplete_peptides
                )
              }

              # accumulate sample NA counts
              if (!is.null(ebm_result_up$incomplete_peptide_nas)) {
                for (sample_name in names(ebm_result_up$incomplete_peptide_nas)) {
                  if (!sample_name %in% names(ebm_stats$sample_na_counts)) {
                    ebm_stats$sample_na_counts[[sample_name]] <- 0
                  }
                  ebm_stats$sample_na_counts[[sample_name]] <-
                    ebm_stats$sample_na_counts[[sample_name]] + ebm_result_up$incomplete_peptide_nas[[sample_name]]
                }
              }
            }

            if (ebm_result_up$status == "success") {
              ebm_stats$n_success <- ebm_stats$n_success + 1
              results_by_domain_ebm[[paste0(domain_id, "_UP")]] <- data.frame(
                DomainID = domain_id,
                DomainName = up_peptides$DomainName[1],
                ProteinName = up_peptides$ProteinName[1],
                Label = contrast,
                Direction = "UP",
                n_peptides = n_up,
                combined_pvalue = ebm_result_up$P_test,
                mean_log2FC = mean(up_peptides$log2FC),
                median_log2FC = median(up_peptides$log2FC),
                method = "EBM",
                stringsAsFactors = FALSE
              )
            } else {
              ebm_stats[[paste0("n_failed_", ebm_result_up$status)]] <-
                ebm_stats[[paste0("n_failed_", ebm_result_up$status)]] + 1
            }
          }

          # Cauchy method
          cauchy_stats$n_attempted <- cauchy_stats$n_attempted + 1
          tryCatch(
            {
              up_pvals_one_tailed <- up_peptides$adj.pvalue / 2
              cauchy_combined_up <- ACAT(up_pvals_one_tailed)
              cauchy_stats$n_success <- cauchy_stats$n_success + 1

              results_by_domain_cauchy[[paste0(domain_id, "_UP")]] <- data.frame(
                DomainID = domain_id,
                DomainName = up_peptides$DomainName[1],
                ProteinName = up_peptides$ProteinName[1],
                Label = contrast,
                Direction = "UP",
                n_peptides = n_up,
                combined_pvalue = cauchy_combined_up,
                mean_log2FC = mean(up_peptides$log2FC),
                median_log2FC = median(up_peptides$log2FC),
                method = "cauchy",
                stringsAsFactors = FALSE
              )
            },
            error = function(e) {
              cauchy_stats$n_failed <- cauchy_stats$n_failed + 1
              if (verbose) {
                cat(sprintf(" > domain %s UP: cauchy error: %s\n", domain_id, e$message))
              }
            }
          )
        }
      }

      # process DOWN direction #
      ##########################
      if (n_down > 0) {
        direction_stats$n_total_instances <- direction_stats$n_total_instances + 1

        # single peptide DOWN
        if (n_down == 1) {
          direction_stats$n_single_peptide <- direction_stats$n_single_peptide + 1

          one_tailed_pval <- down_peptides$adj.pvalue[1] / 2
          single_down_result <- data.frame(
            DomainID = domain_id,
            DomainName = down_peptides$DomainName[1],
            ProteinName = down_peptides$ProteinName[1],
            Label = contrast,
            Direction = "DOWN",
            n_peptides = 1,
            combined_pvalue = one_tailed_pval,
            mean_log2FC = down_peptides$log2FC[1],
            median_log2FC = down_peptides$log2FC[1],
            method = "single_peptide",
            stringsAsFactors = FALSE
          )
          results_by_domain_ebm[[paste0(domain_id, "_DOWN")]] <- single_down_result
          results_by_domain_cauchy[[paste0(domain_id, "_DOWN")]] <- single_down_result

          # multi-peptide DOWN
        } else if (n_down >= 2) {
          direction_stats$n_multi_peptide <- direction_stats$n_multi_peptide + 1

          # EBM method
          ebm_stats$n_attempted <- ebm_stats$n_attempted + 1
          ebm_result_down <- process_direction_ebm(
            peptides = down_peptides,
            replicate_data = contrast_replicates,
            direction = "DOWN",
            verbose = verbose,
            domain_id = domain_id
          )

          if (!is.null(ebm_result_down$status)) {
            if (ebm_result_down$n_incomplete > 0) {
              ebm_stats$n_peptides_removed <- ebm_stats$n_peptides_removed + ebm_result_down$n_incomplete

              if (!is.null(ebm_result_down$incomplete_peptides)) {
                ebm_stats$unique_incomplete_peptides <- c(
                  ebm_stats$unique_incomplete_peptides,
                  ebm_result_down$incomplete_peptides
                )
              }

              if (!is.null(ebm_result_down$incomplete_peptide_nas)) {
                for (sample_name in names(ebm_result_down$incomplete_peptide_nas)) {
                  if (!sample_name %in% names(ebm_stats$sample_na_counts)) {
                    ebm_stats$sample_na_counts[[sample_name]] <- 0
                  }
                  ebm_stats$sample_na_counts[[sample_name]] <-
                    ebm_stats$sample_na_counts[[sample_name]] + ebm_result_down$incomplete_peptide_nas[[sample_name]]
                }
              }
            }

            if (ebm_result_down$status == "success") {
              ebm_stats$n_success <- ebm_stats$n_success + 1
              results_by_domain_ebm[[paste0(domain_id, "_DOWN")]] <- data.frame(
                DomainID = domain_id,
                DomainName = down_peptides$DomainName[1],
                ProteinName = down_peptides$ProteinName[1],
                Label = contrast,
                Direction = "DOWN",
                n_peptides = n_down,
                combined_pvalue = ebm_result_down$P_test,
                mean_log2FC = mean(down_peptides$log2FC),
                median_log2FC = median(down_peptides$log2FC),
                method = "EBM",
                stringsAsFactors = FALSE
              )
            } else {
              ebm_stats[[paste0("n_failed_", ebm_result_down$status)]] <-
                ebm_stats[[paste0("n_failed_", ebm_result_down$status)]] + 1
            }
          }

          # Cauchy
          cauchy_stats$n_attempted <- cauchy_stats$n_attempted + 1
          tryCatch(
            {
              down_pvals_one_tailed <- down_peptides$adj.pvalue / 2
              cauchy_combined_down <- ACAT(down_pvals_one_tailed)
              cauchy_stats$n_success <- cauchy_stats$n_success + 1

              results_by_domain_cauchy[[paste0(domain_id, "_DOWN")]] <- data.frame(
                DomainID = domain_id,
                DomainName = down_peptides$DomainName[1],
                ProteinName = down_peptides$ProteinName[1],
                Label = contrast,
                Direction = "DOWN",
                n_peptides = n_down,
                combined_pvalue = cauchy_combined_down,
                mean_log2FC = mean(down_peptides$log2FC),
                median_log2FC = median(down_peptides$log2FC),
                method = "cauchy",
                stringsAsFactors = FALSE
              )
            },
            error = function(e) {
              cauchy_stats$n_failed <- cauchy_stats$n_failed + 1
            }
          )
        }
      }
    }

    # contrast summary
    cat("\n\nDomain-Level Breakdown...\n")
    cat(sprintf("%-42s %27d\n", " > single peptide:", domain_stats$n_single_peptide))
    cat(sprintf("%-42s %27d\n", " > multi peptide:", domain_stats$n_multi_peptide))
    cat(sprintf("%-42s %27d\n", "  > only up:", domain_stats$n_up_only))
    cat(sprintf("%-42s %27d\n", "  > only down:", domain_stats$n_down_only))
    cat(sprintf("%-42s %27d\n", "  > mixed:", domain_stats$n_mixed))

    cat("\nDirection-Level Breakdown...\n")
    cat(sprintf("%-42s %27d\n", " > total direction instances:", direction_stats$n_total_instances))
    cat(sprintf("%-42s %27d\n", "  > single-peptide (passthrough):", direction_stats$n_single_peptide))
    cat(sprintf("%-42s %27d\n", "  > multi-peptide (aggregated):", direction_stats$n_multi_peptide))

    cat("\nEBM Breakdown...\n")
    cat(sprintf("%-42s %27d\n", " > tests attempted:", ebm_stats$n_attempted))

    n_failed_ebm <- ebm_stats$n_failed_less_than_two + ebm_stats$n_failed_ebm_error

    cat(sprintf("%-42s %27s\n", " > tests completed:", sprintf(
      "(%5.1f%%) %5d",
      100 * ebm_stats$n_success / ebm_stats$n_attempted, ebm_stats$n_success
    )))

    cat(sprintf(
      "%-42s %27s\n", " > total tests failed:",
      sprintf("(%5.1f%%) %5d", 100 * n_failed_ebm / ebm_stats$n_attempted, n_failed_ebm)
    ))

    cat(sprintf(
      "%-42s %27s\n", "  > tests with less than two peptides",
      sprintf("(%5.1f%%) %5d", 100 * ebm_stats$n_failed_less_than_two / ebm_stats$n_attempted, ebm_stats$n_failed_less_than_two)
    ))

    cat(sprintf(
      "%-42s %27s\n", "  > tests failed (EBM error):",
      sprintf("(%5.1f%%) %5d", 100 * ebm_stats$n_failed_ebm_error / ebm_stats$n_attempted, ebm_stats$n_failed_ebm_error)
    ))

    cat("\nMissing Data Summary (EBM only)...\n")
    n_unique_incomplete <- length(unique(ebm_stats$unique_incomplete_peptides))
    cat(sprintf("%-42s %27d\n", " > peptides with missing data:", n_unique_incomplete))

    if (length(ebm_stats$sample_na_counts) > 0) {
      total_nas <- sum(unlist(ebm_stats$sample_na_counts))
      cat(sprintf("%-42s %27d\n", " > total NAs:", total_nas))
      cat(" > NAs by sample:\n")
      sample_counts <- sort(unlist(ebm_stats$sample_na_counts), decreasing = TRUE)
      for (sample in names(sample_counts)) {
        pct_of_total <- 100 * sample_counts[[sample]] / sum(sample_counts)
        cat(sprintf("   %-39s %27s\n", sample, sprintf("(%5.1f%%) %5d", pct_of_total, sample_counts[[sample]])))
      }
    }

    # cat(sprintf("%-42s %27s\n", "  > failed because less than two peptides:", ebm_stats$n_failed_less_than_two))

    # cauchy summary
    cat("\nCauchy Breakdown...\n")
    cat(sprintf("%-42s %27d\n", " > tests attempted:", cauchy_stats$n_attempted))

    n_failed_cauchy <- cauchy_stats$n_failed

    cat(sprintf("%-42s %27s\n", " > tests completed:", sprintf(
      "(%5.1f%%) %5d",
      100 * cauchy_stats$n_success / cauchy_stats$n_attempted, cauchy_stats$n_success
    )))
    cat(sprintf("%-42s %27s\n\n", " > tests failed:", sprintf(
      "(%5.1f%%) %5d",
      100 * n_failed_cauchy / cauchy_stats$n_attempted, n_failed_cauchy
    )))

    # did everything add up?
    if (direction_stats$n_total_instances !=
      (direction_stats$n_single_peptide + direction_stats$n_multi_peptide)) {
      cat("WARNING: Direction counts don't add up!\n")
    }

    if (ebm_stats$n_attempted != direction_stats$n_multi_peptide) {
      cat("WARNING: EBM tests don't match multi-peptide direction instances!\n")
    }

    # store results by contrast
    if (length(results_by_domain_ebm) > 0) {
      contrast_df_ebm <- bind_rows(results_by_domain_ebm)
      results_by_contrast_ebm[[contrast]] <- contrast_df_ebm

      n_sig_ebm <- sum(contrast_df_ebm$combined_pvalue < pval_cutoff &
        abs(contrast_df_ebm$mean_log2FC) > fc_cutoff, na.rm = TRUE)

      cat(sprintf(
        "%-42s %27s\n", sprintf("EBM Significant (p < %.2f, |FC| > %.1f):", pval_cutoff, fc_cutoff),
        sprintf("(%5.1f%%) %d / %d", 100 * n_sig_ebm / nrow(contrast_df_ebm), n_sig_ebm, nrow(contrast_df_ebm))
      ))
    }

    if (length(results_by_domain_cauchy) > 0) {
      contrast_df_cauchy <- bind_rows(results_by_domain_cauchy)
      results_by_contrast_cauchy[[contrast]] <- contrast_df_cauchy

      n_sig_cauchy <- sum(contrast_df_cauchy$combined_pvalue < pval_cutoff &
        abs(contrast_df_cauchy$mean_log2FC) > fc_cutoff, na.rm = TRUE)

      cat(sprintf(
        "%-42s %27s\n\n", sprintf("Cauchy Significant (p < %.2f, |FC| > %.1f):", pval_cutoff, fc_cutoff),
        sprintf("(%5.1f%%) %d / %d", 100 * n_sig_cauchy / nrow(contrast_df_cauchy), n_sig_cauchy, nrow(contrast_df_cauchy))
      ))
    }

    cat("\n")
    global_domain_stats[[contrast]] <- domain_stats
    global_ebm_stats[[contrast]] <- ebm_stats
    global_cauchy_stats[[contrast]] <- cauchy_stats
  }

  # global summary
  cat("All Contrast Summary\n")
  cat(rep("=", 70), "\n", sep = "")

  # aggregate method stats
  total_ebm_tests <- sum(sapply(global_ebm_stats, function(x) x$n_attempted))
  total_ebm_success <- sum(sapply(global_ebm_stats, function(x) x$n_success))

  total_cauchy_tests <- sum(sapply(global_cauchy_stats, function(x) x$n_attempted))
  total_cauchy_success <- sum(sapply(global_cauchy_stats, function(x) x$n_success))

  cat(sprintf("%-48s %21d\n", "Aggregation count (multi-peptide per direction):", total_ebm_tests))
  cat(sprintf(
    "%-42s %27s\n", " > EBM success rate:",
    sprintf(
      "(%.1f%%) %d / %d", 100 * total_ebm_success / total_ebm_tests,
      total_ebm_success, total_ebm_tests
    )
  ))

  cat(sprintf(
    "%-42s %27s\n\n", " > Cauchy success rate:",
    sprintf(
      "(%.1f%%) %d / %d", 100 * total_cauchy_success / total_cauchy_tests,
      total_cauchy_success, total_cauchy_tests
    )
  ))

  # final results
  final_results_ebm <- bind_rows(results_by_contrast_ebm)
  final_results_cauchy <- bind_rows(results_by_contrast_cauchy)

  if (nrow(final_results_ebm) > 0) {
    final_results_ebm <- final_results_ebm %>%
      arrange(combined_pvalue)      
  }

  if (nrow(final_results_cauchy) > 0) {
    final_results_cauchy <- final_results_cauchy %>%
      arrange(combined_pvalue)
  }

  # return list of results and stats
  return(list(
    ebm = final_results_ebm,
    cauchy = final_results_cauchy,
    stats_ebm = global_ebm_stats,
    stats_cauchy = global_cauchy_stats
  ))
}
