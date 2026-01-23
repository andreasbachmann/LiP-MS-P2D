# groups by identical statistics (same p-value, FC, n_peptides = same peptides)
# keep shortest domain name as representative
deduplicate_results <- function(results) {
  results %>%
    group_by(ProteinName, Label, Direction, mean_log2FC, combined_pvalue, n_peptides) %>%
    summarize(
      DomainName = DomainName[which.min(nchar(DomainName))],
      DomainID = first(DomainID),
      method = first(method),
      .groups = "drop"
    ) %>%
    select(
      DomainID, DomainName, ProteinName, Label, Direction,
      n_peptides, combined_pvalue, mean_log2FC, median_log2FC, method
    ) %>%
    arrange(combined_pvalue) %>%
    as.data.table()
}