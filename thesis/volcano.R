# volcano.R
# creates volcano plot for domain-level results

# install missing packages if needed
if (!require("patchwork", quietly = TRUE)) install.packages("patchwork")
if (!require("ggrepel", quietly = TRUE)) install.packages("ggrepel")

# libraries
library(dplyr)
library(tidyr)
library(data.table)
library(patchwork)
library(ggplot2)
library(ggrepel)


# load domain-level results
cauchy_i <- fread("results/domain_level_results_cauchy_dedup.csv")
ebm_i <- fread("results/domain_level_results_ebm_dedup.csv")

# create one tailed volcano plot
create_onetail_volcano <- function(data,
                                   tail = "right",
                                   contrast = NULL,
                                   title = NULL,
                                   significance_level = 0.01,
                                   log2fc_cutoff = 1,
                                   n_top_labels = 20,
                                   exclude_single = FALSE) {
  data_name <- deparse(substitute(data))

  # pick the right column
  pval_col <- ifelse(tail == "right", "adj_p_right", "adj_p_left")

  data <- data %>%
    filter(!is.na(.data[[pval_col]]))

  if (exclude_single) {
    data <- data %>% filter(method != "single_peptide")
  }

  data <- data %>% mutate(pval = .data[[pval_col]])

  if (!is.null(contrast)) {
    data <- data %>% filter(Label == contrast)
  }

  sig_data <- data %>%
    filter(abs(mean_log2FC) > log2fc_cutoff, pval < significance_level) %>%
    arrange(pval)

  # count how many are on the "expected" vs "unexpected" side
  if (tail == "right") {
    n_expected <- sum(sig_data$mean_log2FC > 0)
    n_unexpected <- sum(sig_data$mean_log2FC <= 0)
    tail_label <- "Right tail (UP)"
  } else {
    n_expected <- sum(sig_data$mean_log2FC < 0)
    n_unexpected <- sum(sig_data$mean_log2FC >= 0)
    tail_label <- "Left tail (DOWN)"
  }

  if (is.null(title)) {
    method_label <- case_when(
      grepl("ebm", data_name, ignore.case = TRUE) ~ "EBM",
      grepl("cauchy", data_name, ignore.case = TRUE) ~ "CCT",
      TRUE ~ data_name
    )
    contrast_label <- ifelse(contrast == "Asernate vs Normal", "Arsenate vs Normal", contrast)
    title <- sprintf("%s - %s - %s", method_label, contrast_label, tail_label)
  }

  label_data <- sig_data %>%
    arrange(pval) %>%
    head(n_top_labels)

  p <- ggplot(data, aes(x = mean_log2FC, y = -log10(pval))) +
    geom_point(aes(color = case_when(
      pval < significance_level & abs(mean_log2FC) >= log2fc_cutoff ~ "p-value & log2FC",
      pval < significance_level & abs(mean_log2FC) < log2fc_cutoff ~ "p-value only",
      pval >= significance_level & abs(mean_log2FC) >= log2fc_cutoff ~ "log2FC only",
      TRUE ~ "NS"
    )), alpha = 0.3) +
    scale_color_manual(
      values = c(
        "NS" = "grey",
        "p-value only" = "#7570b3",
        "log2FC only" = "#1b9e77",
        "p-value & log2FC" = "#d95f02"
      ),
      breaks = c("NS", "p-value only", "log2FC only", "p-value & log2FC")
    ) +
    geom_hline(yintercept = -log10(significance_level), linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed", alpha = 0.5) +
    geom_text_repel(
      data = label_data,
      aes(label = DomainName),
      segment.size = 0.2,
      size = 3,
      max.overlaps = 30,
      min.segment.length = 0
    ) +
    theme_bw() +
    labs(
      x = "mean(log2FC)",
      y = sprintf("-log10(%s)", pval_col),
      color = "Significance",
      title = title,
      subtitle = sprintf(
        "Significant: %d (of %d total domains)",
        nrow(sig_data),
        nrow(data)
      )
    ) +
    theme(legend.position = "top")

  return(p)
}

# create full volcano plot
create_volcano_plot <- function(data,
                                contrast = NULL,
                                title = NULL,
                                significance_level = 0.01,
                                log2fc_cutoff = 1,
                                n_top_pval = 10,
                                n_top_fc = 10,
                                exclude_single = FALSE) {
  data_name <- deparse(substitute(data))

  data <- data %>%
    filter(!is.na(adj_p_right), !is.na(adj_p_left))

  if (exclude_single) {
    data <- data %>% filter(method != "single_peptide")
  }

  data <- data %>%
    mutate(
      Direction = ifelse(mean_log2FC >= 0, "UP", "DOWN"),
      combined_pvalue = ifelse(mean_log2FC >= 0, adj_p_right, adj_p_left)
    )

  if (!is.null(contrast)) {
    data <- data %>% filter(Label == contrast)
  }

  sig_data <- data %>%
    filter(abs(mean_log2FC) > log2fc_cutoff, combined_pvalue < significance_level) %>%
    arrange(combined_pvalue)

  if (is.null(title)) {
    method_label <- case_when(
      grepl("ebm", data_name, ignore.case = TRUE) ~ "EBM",
      grepl("cauchy", data_name, ignore.case = TRUE) ~ "CCT",
      TRUE ~ data_name
    )
    contrast_label <- ifelse(contrast == "Asernate vs Normal", "Arsenate vs Normal", contrast)
    title <- sprintf("%s - %s", method_label, contrast_label)
  }

  # label top n by pvalue
  top_pval <- sig_data %>%
    arrange(combined_pvalue) %>%
    head(n_top_pval)

  # label top n by abs(log2fc)
  top_fc <- sig_data %>%
    arrange(desc(abs(mean_log2FC))) %>%
    head(n_top_fc)

  label_data <- bind_rows(top_pval, top_fc) %>%
    distinct(DomainName, .keep_all = TRUE)

  p <- ggplot(data, aes(x = mean_log2FC, y = -log10(combined_pvalue))) +
    geom_point(aes(color = case_when(
      combined_pvalue < significance_level & abs(mean_log2FC) >= log2fc_cutoff ~ "p-value & log2FC",
      combined_pvalue < significance_level & abs(mean_log2FC) < log2fc_cutoff ~ "p-value only",
      combined_pvalue >= significance_level & abs(mean_log2FC) >= log2fc_cutoff ~ "log2FC only",
      TRUE ~ "NS"
    )), alpha = 0.3) +
    scale_color_manual(
      values = c(
        "NS" = "grey",
        "p-value only" = "#7570b3",
        "log2FC only" = "#1b9e77",
        "p-value & log2FC" = "#d95f02"
      ),
      breaks = c("NS", "p-value only", "log2FC only", "p-value & log2FC")
    ) +
    geom_hline(yintercept = -log10(significance_level), linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed", alpha = 0.5) +
    geom_text_repel(
      data = label_data,
      aes(label = DomainName),
      segment.size = 0.2,
      size = 2.5,
      max.overlaps = 20,
      box.padding = 0.5,
      point.padding = 0.2,
      min.segment.length = 0
    ) +
    theme_bw() +
    labs(
      x = "mean(log2FC)",
      y = "-log10(p-value)",
      color = "Significance",
      title = title,
      subtitle = sprintf(
        "Significant domains: %d (down: %d, up: %d, of %d total)",
        nrow(sig_data),
        sum(sig_data$Direction == "DOWN"),
        sum(sig_data$Direction == "UP"),
        nrow(data)
      )
    ) +
    theme(legend.position = "top")

  return(p)
}


# FIGURE 4: Volcano CCT Arsenate vs Normal
create_volcano_plot(cauchy_i, contrast = "Asernate vs Normal", exclude_single = TRUE)

# FIGURE 5: Volcano EBM Arsenate vs Normal
create_volcano_plot(ebm_i, contrast = "Asernate vs Normal", exclude_single = TRUE)

# FIGURE 6: Volcano CCT Heat vs Normal
create_volcano_plot(cauchy_i, contrast = "Heat vs Normal", exclude_single = TRUE)