# groups by identical statistics (same p-value, log2FC, n_peptides, ... = same peptides)
# keep shortest domain name as representative
deduplicate_results <- function(results, global_mode = FALSE) {
  if (global_mode) {
    # global mode
    deduped <- results %>%
      group_by(Label, mean_log2FC, median_log2FC, adj_p_right, adj_p_left, n_peptides) %>%
      slice_min(nchar(DomainName), n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      select(
        DomainName, ProteinName, Label, n_peptides,
        adj_p_right, adj_p_left, mean_log2FC, median_log2FC, method,
        n_proteins, n_domain_instances
      )
  } else {
    # instance mode
    deduped <- results %>%
      group_by(ProteinName, Label, mean_log2FC, median_log2FC, adj_p_right, adj_p_left, n_peptides) %>%
      slice_min(nchar(DomainName), n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      select(
        DomainID, DomainName, ProteinName, Label, n_peptides,
        adj_p_right, adj_p_left, mean_log2FC, median_log2FC, method
      )
  }

  deduped %>%
    arrange(pmin(adj_p_right, adj_p_left)) %>%
    as.data.table()
}