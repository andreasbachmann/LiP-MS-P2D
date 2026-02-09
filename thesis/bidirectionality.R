# bidirectionality.R
# analyze bidirectionality of domain significance

# libs
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

rm(list = ls())

# load data
ebm_i <- fread("results/domain_level_results_ebm_dedup.csv")
cauchy_i <- fread("results/domain_level_results_cauchy_dedup.csv")


# bidirectionality analysis
calc_bidirectionality <- function(df, label, p_threshold = 0.05,
                                  mode = "instance", exclude_single = FALSE) {
  # set grouping column
  id_col <- ifelse(mode == "global", "DomainName", "DomainID")

  df_label <- df %>%
    filter(Label == label, !is.na(adj_p_right), !is.na(adj_p_left))

  if (exclude_single) {
    df_label <- df_label %>% filter(method != "single_peptide")
  }

  n_total <- n_distinct(df_label[[id_col]])

  bidirectional <- df_label %>%
    filter(adj_p_right < p_threshold, adj_p_left < p_threshold)

  n_bidir <- nrow(bidirectional)

  tibble(
    total_domains = n_total,
    n_bidirectional = n_bidir,
    fraction = n_bidir / n_total
  )
}

get_top_bidirectional <- function(df, label, n = 5, p_threshold = 0.05,
                                  mode = "instance", exclude_single = FALSE) {
  id_col <- ifelse(mode == "global", "DomainName", "DomainID")

  df_label <- df %>%
    filter(Label == label, !is.na(adj_p_right), !is.na(adj_p_left))

  if (exclude_single) {
    df_label <- df_label %>% filter(method != "single_peptide")
  }

  df_label %>%
    filter(adj_p_right < p_threshold, adj_p_left < p_threshold) %>%
    mutate(mean_adj_p = (adj_p_right + adj_p_left) / 2) %>%
    arrange(mean_adj_p) %>%
    head(n)
}


contrast <- "Asernate vs Normal"

# fractions
calc_bidirectionality(ebm_i, contrast, mode = "instance")
calc_bidirectionality(cauchy_i, contrast, mode = "instance")
calc_bidirectionality(ebm_i, contrast, mode = "instance", exclude_single = TRUE)
calc_bidirectionality(cauchy_i, contrast, mode = "instance", exclude_single = TRUE)

# top hits
get_top_bidirectional(ebm_i, contrast, n = 5, mode = "instance")
get_top_bidirectional(cauchy_i, contrast, n = 5, mode = "instance")

# bar plot
plot_bidirectionality <- function(ebm_df, cauchy_df, label,
                                  mode = "instance", exclude_single = FALSE) {
  subtitle_parts <- c(
    gsub("Asernate", "Arsenate", label),
    mode,
    if (exclude_single) "no single-peptide" else NULL
  )

  plot_data <- tibble(
    method = c("EBM", "CCT"),
    fraction_pct = c(
      calc_bidirectionality(ebm_df, label, mode = mode, exclude_single = exclude_single)$fraction,
      calc_bidirectionality(cauchy_df, label, mode = mode, exclude_single = exclude_single)$fraction
    ) * 100
  )

  ggplot(plot_data, aes(x = method)) +
    geom_col(aes(y = 100), fill = "#b9b9b9", width = 0.6) +
    geom_col(aes(y = fraction_pct, fill = method), width = 0.6) +
    geom_text(aes(y = fraction_pct + 3, label = sprintf("%.1f%%", fraction_pct)), size = 4) +
    scale_fill_manual(values = c("EBM" = "#4DAF4A", "CCT" = "#377EB8")) +
    labs(
      title = "Domains with Bidirectional Significance",
      subtitle = paste(subtitle_parts, collapse = ", "),
      x = NULL, y = "Fraction of Domains (%)"
    ) +
    ylim(0, 105) +
    theme_minimal() +
    theme(legend.position = "none")
}

plot_bidirectionality(ebm_i, cauchy_i, contrast, mode = "instance")
plot_bidirectionality(ebm_i, cauchy_i, contrast, mode = "instance", exclude_single = TRUE)