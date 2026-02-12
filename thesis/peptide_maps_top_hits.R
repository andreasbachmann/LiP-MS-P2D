# peptide_maps_top_hits.R
# check top hits and generate peptide maps

# install missing packages if needed
if (!require("patchwork", quietly = TRUE)) install.packages("patchwork")
if (!require("ggrepel", quietly = TRUE)) install.packages("ggrepel")

# load libraries
library(dplyr)
library(ggplot2)
library(data.table)
library(tidyverse)
library(patchwork)
library(ggrepel)


# load data
df_binding <- fread("data/domain_annotations_consolidated.csv")
df_peptides <- fread("results/peptide_domain_mapping.csv")
cauchy_i <- fread("results/domain_level_results_cauchy_dedup.csv")
ebm_i <- fread("results/domain_level_results_ebm_dedup.csv")

# get pvalue for correct direction
add_combined <- function(df) {
  df %>%
    filter(!is.na(adj_p_right), !is.na(adj_p_left)) %>%
    mutate(
      Direction = ifelse(mean_log2FC >= 0, "UP", "DOWN"),
      combined_pvalue = ifelse(mean_log2FC >= 0, adj_p_right, adj_p_left)
    )
}

cauchy_i <- add_combined(cauchy_i)
ebm_i <- add_combined(ebm_i)


# helper functions
get_hits <- function(df, ids, n = 3, top = TRUE) {
  df %>%
    filter(Label == "Asernate vs Normal", DomainID %in% ids) %>%
    arrange(if (top) combined_pvalue else desc(combined_pvalue)) %>%
    head(n)
}


# plot: top 20 bar chart
plot_top20 <- function(df, label, subtitle_method, id_col = "DomainName", max_peptides = NULL) {
  top <- df %>%
    filter(Label == label) %>%
    arrange(combined_pvalue) %>%
    slice_head(n = 20) %>%
    arrange(mean_log2FC) %>%
    mutate(
      plot_label = factor(.data[[id_col]], levels = .data[[id_col]]),
      n_peptides_display = if (!is.null(max_peptides)) pmin(n_peptides, max_peptides) else n_peptides,
      capped = if (!is.null(max_peptides)) n_peptides > max_peptides else FALSE
    )

  p <- ggplot(top, aes(y = plot_label, x = n_peptides_display, fill = mean_log2FC)) +
    geom_col(width = 0.7, color = "grey70", linewidth = 0.3) +
    scale_fill_gradient2(
      low = "#2673d8", mid = "white", high = "#d73027",
      midpoint = 0, name = "mean log2FC"
    ) +
    theme_bw() +
    labs(
      x = "Number of peptides", y = NULL,
      title = "Top 20 significant domains",
      subtitle = sprintf(
        "%s - %s", subtitle_method,
        ifelse(label == "Asernate vs Normal", "Arsenate vs Normal", label)
      )
    )

  # label capped bars with actual count
  if (!is.null(max_peptides) && any(top$capped)) {
    p <- p +
      geom_text(
        data = top %>% filter(capped),
        aes(label = n_peptides),
        hjust = -0.2, size = 3
      ) +
      xlim(0, max_peptides * 1.15)
  }

  p
}


# plot: protein peptide map
plot_protein_peptide_map <- function(
  protein_id, df_peptides, df_annotations,
  contrast = "Asernate vs Normal", domain_name = NULL,
  adj_pval_threshold = 0.05, domain_height = 0.2, domain_gap = 0.1,
  save_plot = FALSE, output_dir = "~/masterProject/figures/final/"
) {
  df_protein <- df_peptides %>%
    filter(
      ProteinName == protein_id, Label == contrast,
      !is.na(log2FC), is.finite(log2FC),
      !is.na(p_right), !is.na(p_left)
    ) %>%
    mutate(
      directional_pval = ifelse(log2FC >= 0, p_right, p_left),
      significant = directional_pval < adj_pval_threshold
    )

  if (nrow(df_protein) == 0) {
    warning(paste("No peptides for", protein_id))
    return(NULL)
  }

  ymin_pep <- min(df_protein$log2FC)

  if (!is.null(domain_name)) {
    domains_target <- df_annotations %>%
      filter(ProteinName == protein_id, DomainName %in% domain_name) %>%
      select(DomainName, DomainStart, DomainEnd) %>%
      distinct() %>%
      arrange(DomainStart) %>%
      mutate(
        idx = row_number(),
        ymax_dom = ymin_pep - domain_gap - (idx - 1) * (domain_height + domain_gap),
        ymin_dom = ymax_dom - domain_height
      )
  } else {
    domains_target <- data.frame()
  }

  p <- ggplot() +
    geom_segment(
      data = df_protein,
      aes(x = StartPos, xend = EndPos, y = log2FC, yend = log2FC, color = significant),
      linewidth = 1.5, alpha = 0.7
    ) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
    scale_color_manual(
      values = c("FALSE" = "grey50", "TRUE" = "#d95f02"),
      labels = c("NS", "Significant")
    ) +
    labs(
      x = "AA position", y = "log2FC",
      title = protein_id,
      color = ""
    ) +
    theme_bw() +
    theme(legend.position = "top")

  if (nrow(domains_target) > 0) {
    p <- p +
      geom_rect(
        data = domains_target,
        aes(
          xmin = DomainStart, xmax = DomainEnd,
          ymin = ymin_dom, ymax = ymax_dom, fill = DomainName
        ),
        alpha = 0.6
      ) +
      scale_fill_brewer(palette = "Set2", name = "Domain")
  }

  if (save_plot) {
    filename <- paste0(output_dir, "peptide_map_", protein_id, ".png")
    ggsave(filename, plot = p, width = 10, height = 4, dpi = 300)
    message("Saved: ", filename)
  }

  return(p)
}

plot_multi_protein_map <- function(
  protein_ids, df_peptides, df_annotations,
  contrast = "Asernate vs Normal", domain_name = NULL,
  adj_pval_threshold = 0.05, domain_height = 0.2, domain_gap = 0.1,
  ncol = 1, save_plot = FALSE,
  output_dir = "~/masterProject/figures/final/tophits/"
) {
  plots <- lapply(protein_ids, function(pid) {
    plot_protein_peptide_map(
      protein_id = pid, df_peptides = df_peptides,
      df_annotations = df_annotations, contrast = contrast,
      domain_name = domain_name, adj_pval_threshold = adj_pval_threshold,
      domain_height = domain_height, domain_gap = domain_gap
    )
  })
  plots <- plots[!sapply(plots, is.null)]
  if (length(plots) == 0) {
    return(NULL)
  }

  combined <- wrap_plots(plots, ncol = ncol) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  if (save_plot) {
    protein_str <- paste(protein_ids, collapse = "_")
    domain_str <- if (!is.null(domain_name)) {
      paste0("_", gsub("/", "-", gsub(" ", "-", gsub(",", "", paste(domain_name, collapse = "_")))))
    } else {
      ""
    }
    filename <- paste0(output_dir, "peptide_map_", gsub(" ", "_", contrast), "_", protein_str, domain_str, ".png")
    ggsave(filename, plot = combined, width = 12, height = 4 * ceiling(length(plots) / ncol), dpi = 300)
    message("Saved: ", filename)
  }

  return(combined)
}

plot_hits <- function(hits_df, df_peptides, df_binding,
                      title_prefix = "", subtitle = NULL,
                      save_plot = FALSE, ncol = 2) {
  protein_domains <- hits_df %>%
    group_by(ProteinName) %>%
    summarise(domains = list(unique(DomainName)), .groups = "drop")

  plots <- lapply(seq_len(nrow(protein_domains)), function(i) {
    plot_multi_protein_map(
      protein_ids = protein_domains$ProteinName[i],
      df_peptides = df_peptides, df_annotations = df_binding,
      domain_name = protein_domains$domains[[i]]
    )
  })
  plots <- plots[!sapply(plots, is.null)]
  if (length(plots) == 0) {
    return(NULL)
  }

  combined <- wrap_plots(plots, ncol = ncol) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = if (title_prefix != "") title_prefix else NULL,
      subtitle = subtitle,
      theme = theme(
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 12)
      )
    ) &
    theme(legend.position = "bottom")

  if (save_plot && title_prefix != "") {
    filename <- sprintf(
      "~/masterProject/figures/final/tophits/%s.png",
      gsub(" ", "_", tolower(title_prefix))
    )
    ggsave(filename,
      plot = combined,
      width = ifelse(ncol == 1, 12, 10 * ncol),
      height = ifelse(ncol == 1, 4, 4) * ceiling(length(plots) / ncol),
      dpi = 300
    )
    message("Saved: ", filename)
  }

  combined
}


# generate figures

# FIGURE 7: Top 20 Hits CCT
plot_top20(cauchy_i, "Asernate vs Normal", "CCT (instance)", id_col = "DomainID")
# ggsave("~/masterProject/figures/final/top20_cauchy_instance.png", width = 10, height = 6, dpi = 300)

# FIGURE 8: Top 20 Hits EBM
plot_top20(ebm_i, "Asernate vs Normal", "EBM (instance)", id_col = "DomainID")
# ggsave("~/masterProject/figures/final/top20_ebm_instance.png", width = 10, height = 6, dpi = 300)


# method agreement
sig_ebm <- ebm_i %>%
  filter(Label == "Asernate vs Normal", combined_pvalue < 0.05) %>%
  pull(DomainID) %>%
  unique()

sig_cauchy <- cauchy_i %>%
  filter(Label == "Asernate vs Normal", combined_pvalue < 0.05) %>%
  pull(DomainID) %>%
  unique()

ebm_only <- setdiff(sig_ebm, sig_cauchy)
cauchy_only <- setdiff(sig_cauchy, sig_ebm)
consensus <- intersect(sig_ebm, sig_cauchy)

# method comparison bar chart prep
combined_hits <- bind_rows(
  get_hits(ebm_i, ebm_only, top = TRUE) %>% mutate(Category = "Top 3 EBM-only"),
  get_hits(cauchy_i, cauchy_only, top = TRUE) %>% mutate(Category = "Top 3 CCT-only"),
  get_hits(cauchy_i, consensus, top = TRUE) %>% mutate(Category = "Top 3 Consensus"),
  get_hits(ebm_i, ebm_only, top = FALSE) %>% mutate(Category = "Bottom 3 EBM-only"),
  get_hits(cauchy_i, cauchy_only, top = FALSE) %>% mutate(Category = "Bottom 3 CCT-only"),
  get_hits(cauchy_i, consensus, top = FALSE) %>% mutate(Category = "Bottom 3 Consensus")
) %>%
  mutate(Category = factor(Category, levels = c(
    "Top 3 EBM-only", "Top 3 CCT-only", "Top 3 Consensus",
    "Bottom 3 EBM-only", "Bottom 3 CCT-only", "Bottom 3 Consensus"
  )))

# FIGURE 11
top_hits <- combined_hits %>% filter(grepl("^Top", Category))
ggplot(top_hits, aes(y = reorder(DomainID, n_peptides), x = n_peptides, fill = mean_log2FC)) +
  geom_col(width = 0.7, color = "grey70", linewidth = 0.3) +
  facet_wrap(~Category, scales = "free", ncol = 1) +
  scale_fill_gradient2(low = "#2673d8", mid = "white", high = "#d73027", midpoint = 0, name = "mean log2FC") +
  theme_bw() +
  labs(
    x = "Number of Peptides", y = NULL,
    title = "Method Comparison: Top Hits",
    subtitle = "Arsenate vs Normal"
  ) +
  theme(legend.position = "right", strip.text = element_text(face = "bold"))
# ggsave("~/masterProject/figures/final/method_comparison_top.png", width = 10, height = 8, dpi = 300)

# FIGURE 12
bottom_hits <- combined_hits %>% filter(grepl("^Bottom", Category))
ggplot(bottom_hits, aes(y = reorder(DomainID, n_peptides), x = n_peptides, fill = mean_log2FC)) +
  geom_col(width = 0.7, color = "grey70", linewidth = 0.3) +
  facet_wrap(~Category, scales = "free", ncol = 1) +
  scale_fill_gradient2(low = "#2673d8", mid = "white", high = "#d73027", midpoint = 0, name = "mean log2FC") +
  theme_bw() +
  labs(
    x = "Number of Peptides", y = NULL,
    title = "Method Comparison: Borderline Hits",
    subtitle = "Arsenate vs Normal"
  ) +
  theme(legend.position = "right", strip.text = element_text(face = "bold"))
# ggsave("~/masterProject/figures/final/method_comparison_bottom.png", width = 10, height = 8, dpi = 300)




# Figure 22: Bottom 3 EBM-only
plot_hits(get_hits(ebm_i, ebm_only, n = 3, top = FALSE),
  df_peptides, df_binding, title_prefix = "Bottom 3 EBM-only",
  subtitle = "Arsenate vs Normal",
  save_plot = TRUE, ncol = 1
)

# Figure 23: Bottom 3 CCT-only
plot_hits(get_hits(cauchy_i, cauchy_only, n = 3, top = TRUE),
  df_peptides, df_binding,
  title_prefix = "Top 3 CCT-only",
  subtitle = "Arsenate vs Normal",
  save_plot = TRUE, ncol = 1
)

# FIGURE 24: Bottom 3 Consensus
plot_hits(get_hits(cauchy_i, consensus, n = 3, top = FALSE),
  df_peptides, df_binding,
  title_prefix = "Bottom 3 Consensus",
  subtitle = "Arsenate vs Normal",
  save_plot = TRUE, ncol = 1
)


# FIGURE 18: Top 5 Bidirectional
bidir_hits <- ebm_i %>%
  filter(Label == "Asernate vs Normal", adj_p_right < 0.05, adj_p_left < 0.05) %>%
  mutate(mean_adj_p = (adj_p_right + adj_p_left) / 2) %>%
  arrange(mean_adj_p) %>%
  head(5)
plot_hits(bidir_hits, df_peptides, df_binding,
  title_prefix = "Top 5 Bidirectionally Significant Domains",
  subtitle = "Arsenate vs Normal (EBM, instance-level)",
  save_plot = TRUE, ncol = 1
)

# FIGURE 15: Top RNA-binding domain CCT
plot_hits(
  cauchy_i %>%
    filter(DomainID == "P08621_RNA recognition motif domain_103_181", Label == "Asernate vs Normal"),
  df_peptides, df_binding,
  title_prefix = "Top-Ranked RNA-Binding Domain Instance (CCT))",
  subtitle = "Arsenate vs Normal",
  save_plot = TRUE, ncol = 1
)

# FIGURE 16: Top RNA-binding domain EBM
plot_hits(
  ebm_i %>%
    filter(DomainID == "P19338_RNA recognition motif domain_4_572_647", Label == "Asernate vs Normal"),
  df_peptides, df_binding,
  title_prefix = "Top-Ranked RNA-Binding Domain Instance (EBM)",
  subtitle = "Arsenate vs Normal",
  save_plot = TRUE, ncol = 1
)

# FIGURE 9: Top ranked domain CCT
plot_hits(
  cauchy_i %>%
    filter(DomainID == "Q15020_LSM-interacting domain_944_962", Label == "Asernate vs Normal"),
  df_peptides, df_binding,
  title_prefix = "Top-Ranked Domain Instance",
  subtitle = "Arsenate vs Normal",
  save_plot = TRUE, ncol = 1
)

# FIGURE 10: annotation redundancy
plot_hits(
  cauchy_i %>%
    filter(
      DomainID %in% c(
        "P28838_Peptidase M17, leucyl aminopeptidase, N-terminal_37_169",
        "P28838_Macro domain-like_33_197",
        "Q5JTZ9_Alanine--tRNA ligase_34_981",
        "Q5JTZ9_Alanyl-tRNA synthetase, class IIc, N-terminal_42_623",
        "Q5JTZ9_Alanine-tRNA ligase, class IIc, anti-codon-binding domain superfamily_285_493"
      ),
      Label == "Asernate vs Normal"
    ),
  df_peptides, df_binding,
  title_prefix = "Annotation Redundancy",
  subtitle = "Arsenate vs Normal",
  save_plot = TRUE, ncol = 1
)


# decent example
plot_hits(
  cauchy_i %>%
    filter(DomainID == "O75822_Eukaryotic translation initiation factor 3 subunit J_12_258", Label == "Asernate vs Normal"),
  df_peptides, df_binding,
  title_prefix = "EXAMPLE",
  subtitle = "Arsenate vs Normal",
  save_plot = FALSE, ncol = 1
)
