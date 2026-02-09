library(dplyr)
library(data.table)
library(ggplot2)
library(scales)

cauchy_i <- fread("results/domain_level_results_cauchy_dedup.csv")
ebm_i <- fread("results/domain_level_results_ebm_dedup.csv")
unique_rbd <- fread("data/unique_RBDs.tsv", sep = "\t")

run_fisher <- function(df, label, direction = NULL, p_cutoff = 0.05) {
  universe <- df %>%
    filter(Label == label, !is.na(adj_p_right), !is.na(adj_p_left)) %>%
    mutate(is_RNA = DomainName %in% unique_rbd$DomainName)

  if (is.null(direction)) {
    sig <- universe %>% filter(adj_p_right < p_cutoff | adj_p_left < p_cutoff)
  } else if (direction == "DOWN") {
    sig <- universe %>% filter(adj_p_left < p_cutoff)
  } else {
    sig <- universe %>% filter(adj_p_right < p_cutoff)
  }

  tbl <- matrix(
    c(
      sum(sig$is_RNA), sum(!sig$is_RNA),
      sum(universe$is_RNA) - sum(sig$is_RNA),
      sum(!universe$is_RNA) - sum(!sig$is_RNA)
    ),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("sig", "non-sig"), c("RNA", "non-RNA"))
  )
  fisher.test(tbl)
}

make_plot_data <- function(df, label, direction = NULL, p_cutoff = 0.05) {
  universe <- df %>%
    filter(Label == label, !is.na(adj_p_right), !is.na(adj_p_left)) %>%
    mutate(is_RNA = DomainName %in% unique_rbd$DomainName)

  if (is.null(direction)) {
    sig <- universe %>% filter(adj_p_right < p_cutoff | adj_p_left < p_cutoff)
    name <- "UP or DOWN"
  } else if (direction == "DOWN") {
    sig <- universe %>% filter(adj_p_left < p_cutoff)
    name <- "DOWN"
  } else {
    sig <- universe %>% filter(adj_p_right < p_cutoff)
    name <- "UP"
  }

  ft <- run_fisher(df, label, direction, p_cutoff)

  tibble(
    Set = c("Universe", "Significant"),
    Fraction_RNA = c(
      sum(universe$is_RNA) / nrow(universe),
      sum(sig$is_RNA) / nrow(sig)
    ),
    Direction = name,
    p_value = ft$p.value,
    OR = ft$estimate
  )
}

# run fisher tests
run_fisher(cauchy_i, "Asernate vs Normal", "DOWN")
run_fisher(cauchy_i, "Asernate vs Normal", "UP")
run_fisher(cauchy_i, "Asernate vs Normal")

run_fisher(ebm_i, "Asernate vs Normal", "DOWN")
run_fisher(ebm_i, "Asernate vs Normal", "UP")
run_fisher(ebm_i, "Asernate vs Normal")

# build plot data
plot_df_cauchy_i <- bind_rows(
  make_plot_data(cauchy_i, "Asernate vs Normal", "DOWN"),
  make_plot_data(cauchy_i, "Asernate vs Normal", "UP"),
  make_plot_data(cauchy_i, "Asernate vs Normal")
)

plot_df_ebm_i <- bind_rows(
  make_plot_data(ebm_i, "Asernate vs Normal", "DOWN"),
  make_plot_data(ebm_i, "Asernate vs Normal", "UP"),
  make_plot_data(ebm_i, "Asernate vs Normal")
)

# plot function to avoid repeating
plot_fisher <- function(plot_df, title_method) {
  ggplot(plot_df, aes(x = Set, y = Fraction_RNA, fill = Set)) +
    geom_col(width = 0.6) +
    facet_wrap(~Direction) +
    scale_y_continuous(labels = percent_format(), limits = c(0, NA)) +
    scale_fill_manual(values = c("Universe" = "#BDBDBD", "Significant" = "#3182BD")) +
    geom_text(
      data = plot_df %>% distinct(Direction, p_value, OR),
      aes(
        x = 1.5,
        y = max(plot_df$Fraction_RNA) * 1.05,
        label = paste0(
          "p = ", ifelse(p_value < 0.001, formatC(p_value, format = "e", digits = 2), round(p_value, 4)),
          "\nOR = ", round(OR, 2)
        )
      ),
      inherit.aes = FALSE, size = 3.5
    ) +
    labs(
      y = "Fraction of RNA-binding domains",
      x = NULL,
      title = "RNA-binding domain enrichment (instance-level)",
      subtitle = paste0(title_method, " (Arsenate vs Normal)")
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none", strip.text = element_text(face = "bold"))
}

# FIGURE 13: RNA enrichment CCT
plot_fisher(plot_df_cauchy_i, "Cauchy Combination Test")
# ggsave("~/masterProject/figures/final/rna_fisher_cauchy.png", width = 6, height = 6)

# FIGURE 14: RNA enrichment EBM
plot_fisher(plot_df_ebm_i, "Empirical Brown's Method")
# ggsave("~/masterProject/figures/final/rna_fisher_ebm.png", width = 6, height = 6)