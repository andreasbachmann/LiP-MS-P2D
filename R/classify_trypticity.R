# function to classify peptides as fully or semi tryptic
classify_trypticity <- function(peptide_results) {
  peptide_results %>%
    mutate(
      trypticity = case_when(
        (NSEMI_TRI & !NTERMINUS) | (CSEMI_TRI & !CTERMINUS) ~ "semi_tryptic",
        TRUE ~ "fully_tryptic"
      )
    )
}