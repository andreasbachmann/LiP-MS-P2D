# determine which peptides belong to which domains based on coordinates
map_peptides_to_domains <- function(peptide_results, domain_annotations, overlap_cutoff = 0.5) {
  
  # extract unique peptides with position
  peptides <- peptide_results %>%
    select(ProteinName, PeptideSequence, FULL_PEPTIDE, StartPos, EndPos) %>%
    distinct()
  
  # perform mapping
  peptide_domains <- peptides %>%
    # join mapping dataframe and pair every peptide with every domain on its protein
    inner_join(domain_annotations, by = "ProteinName", relationship = "many-to-many") %>%
    # filter peptides that dont overlap; must overlap with domain
    filter(StartPos <= DomainEnd & EndPos >= DomainStart) %>%
    # calculate overlapping region fraction
    mutate(
      overlap_start = pmax(StartPos, DomainStart),
      overlap_end = pmin(EndPos, DomainEnd),
      overlap_length = overlap_end - overlap_start + 1,
      peptide_length = EndPos - StartPos + 1,
      overlap_fraction = overlap_length / peptide_length
    ) %>%
    # require percentage overlap?
    filter(overlap_fraction >= overlap_cutoff) %>%
    select(ProteinName, FULL_PEPTIDE, PeptideSequence, StartPos, EndPos, DomainID, DomainName, overlap_fraction)

  return(peptide_domains)
}