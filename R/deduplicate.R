# groups by identical statistics (same p-value, log2FC, n_peptides = same peptides)
# keep shortest domain name as representative
deduplicate_results <- function(results, global_mode = FALSE) {
  if (global_mode) {
    # global mode
    deduped <- results %>%
      group_by(DomainID, Label, Direction, mean_log2FC, median_log2FC, combined_pvalue, n_peptides) %>%
      summarize(
        DomainName = DomainName[which.min(nchar(DomainName))],
        DomainID = first(DomainID),
        ProteinName = first(ProteinName),
        method = first(method),
        n_proteins = first(n_proteins),
        n_domain_instances = first(n_domain_instances),
        .groups = "drop"
      ) %>%
      select(
        DomainName, ProteinName, Label, Direction, n_peptides,
        combined_pvalue, mean_log2FC, median_log2FC, method,
        n_proteins, n_domain_instances
      )
  } else {
    # instance mode
    deduped <- results %>%
      group_by(ProteinName, Label, Direction, mean_log2FC, median_log2FC, combined_pvalue, n_peptides) %>%
      summarize(
        DomainName = DomainName[which.min(nchar(DomainName))],
        DomainID = first(DomainID),
        method = first(method),
        .groups = "drop"
      ) %>%
      select(
        DomainID, DomainName, ProteinName, Label, Direction, n_peptides,
        combined_pvalue, mean_log2FC, median_log2FC, method
      )
  }

  deduped %>%
    arrange(combined_pvalue) %>%
    as.data.table()
}