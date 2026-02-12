# method_agreement.R
# analyze agreement between CCT and EBM

# install missing packages if needed
if (!require("viridis", quietly = TRUE)) install.packages("viridis")

# libraries
library(ggplot2)
library(dplyr)
library(data.table)
library(viridis)

# load instance-level results
ebm_i <- fread("results/domain_level_results_ebm_dedup.csv")
cauchy_i <- fread("results/domain_level_results_cauchy_dedup.csv")

# filter single-peptide, derive combined pvalue
prep <- function(df) {
  df %>%
    filter(method != "single_peptide", !is.na(adj_p_right), !is.na(adj_p_left)) %>%
    mutate(combined_pvalue = ifelse(mean_log2FC >= 0, adj_p_right, adj_p_left))
}

ebm_prep <- prep(ebm_i)
cauchy_prep <- prep(cauchy_i)

# merge
comparison <- inner_join(
  ebm_prep %>% select(DomainID, Label, DomainName, ProteinName,
    n_peptides, mean_log2FC,
    pvalue_ebm = combined_pvalue
  ),
  cauchy_prep %>% select(DomainID, Label, pvalue_cauchy = combined_pvalue),
  by = c("DomainID", "Label")
)

# peptide categories + "All"
comparison <- comparison %>%
  mutate(
    peptide_category = case_when(
      n_peptides == 2 ~ "2 peptides",
      n_peptides <= 5 ~ "3 - 5 peptides",
      n_peptides <= 10 ~ "6 - 10 peptides",
      n_peptides <= 20 ~ "11 - 20 peptides",
      TRUE ~ "> 20 peptides"
    )
  )

comparison_plot <- bind_rows(
  comparison,
  comparison %>% mutate(peptide_category = "All")
) %>%
  mutate(
    peptide_category = factor(peptide_category, levels = c(
      "2 peptides", "3 - 5 peptides", "6 - 10 peptides",
      "11 - 20 peptides", "> 20 peptides", "All"
    ))
  )

# correlation stats per facet
cor_stats <- comparison_plot %>%
  group_by(peptide_category) %>%
  summarise(
    n = n(),
    r = cor(-log10(pvalue_cauchy), -log10(pvalue_ebm), method = "pearson"),
    rho = cor(-log10(pvalue_cauchy), -log10(pvalue_ebm), method = "spearman"),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("n = %s\nr = %.3f\nρ = %.3f", format(n, big.mark = ","), r, rho))

# FIGURE 3: Method Agreement
ggplot(comparison_plot, aes(x = -log10(pvalue_cauchy), y = -log10(pvalue_ebm))) +
  geom_hex(bins = 60) +
  scale_fill_viridis(trans = "log10", name = "Count", option = "plasma") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_text(
    data = cor_stats,
    aes(x = Inf, y = Inf, label = label),
    hjust = 1.1, vjust = 1.5, size = 4, inherit.aes = FALSE
  ) +
  facet_wrap(~peptide_category, ncol = 3) +
  labs(
    title = "CCT vs EBM correlation by peptide count (instance-level)",
    x = "-log10(CCT p-value)",
    y = "-log10(EBM p-value)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "right",
    aspect.ratio = 1
  )
# ggsave("~/masterProject/figures/final/method_agreement.png",width = 10, height = 8, dpi = 300)
