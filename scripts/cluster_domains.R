#!/usr/bin/env Rscript

# load libraries
library(dplyr)
library(data.table)

# parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript cluster_domains.R <input_file> <output_file>")
}

input_file <- args[1]
output_file <- args[2]

# consolidate overlapping domains
cluster_domains <- function(domain_annotations_uncollapsed, max_gap = 5) {
  
  # find different clusters of the same domain on the same protein; to collapse others
  find_domain_clusters <- function(df) {
    
    # only 1 cluster?
    if (nrow(df) == 1) {
      df$ClusterID <- 1
      return(df)
    }
    
    # sort by start position
    df <- df %>% arrange(DomainStart)
    
    # initialize clusters
    clusters <- numeric(nrow(df))
    clusters[1] <- 1
    current_cluster <- 1
    current_end <- df$DomainEnd[1]
    
    # iterate through domains
    for (i in 2:nrow(df)) {
      
      # overlaps with current cluster?
      if (df$DomainStart[i] <= current_end + max_gap) {
        clusters[i] <- current_cluster
        current_end <- max(current_end, df$DomainEnd[i])
      }
      # no overlap? start new cluster
      else {
        current_cluster <- current_cluster + 1
        clusters[i] <- current_cluster
        current_end <- df$DomainEnd[i]
      }
    }
    # write clusterIDs to dataframe and return
    df$ClusterID <- clusters
    return(df)
  }
  
  consolidated <- domain_annotations_uncollapsed %>%
    # group by protein and domain
    group_by(ProteinName, DomainName) %>%
    # apply find_domain to each group
    group_modify(~find_domain_clusters(.x)) %>%
    # group by protein, domain and clusterID
    group_by(ProteinName, DomainName, ClusterID) %>%
    summarize(
      # merge domains in same cluster into one domain
      DomainStart = min(DomainStart),
      DomainEnd = max(DomainEnd),
      # combine IPR_IDs and databases and count merge partners
      IPR_ID = paste(unique(IPR_ID), collapse = ";"),
      Databases = paste(unique(Database), collapse = ";"),
      MergePartners = n(),
      .groups = "drop") %>%
    # create names for multiple instances
    group_by(ProteinName, DomainName) %>%
    arrange(DomainStart) %>%
    mutate(
      n_instances = n(),
      instance = row_number(),
      DomainNameFull = case_when(
        n_instances == 1 ~ DomainName,
        n_instances == 2 ~ paste0(DomainName, "_", c("N-term", "C-term")[instance]),
        TRUE ~ paste0(DomainName, "_", instance)),
      DomainID = paste(ProteinName, DomainNameFull, DomainStart, DomainEnd, sep = "_")) %>%
    ungroup() %>%
    select(ProteinName, IPR_ID, DomainName, Databases, DomainStart, 
           DomainEnd, DomainNameFull, DomainID, ClusterID, MergePartners) %>%
    arrange(ProteinName, DomainName, DomainStart)
  
  return(consolidated)
}

# read input, process, write output
cat("Reading", input_file, "\n")
domains <- fread(input_file)
domains <- domains %>%
  rename(
    ProteinName = V1,
    IPR_ID      = V2,
    DomainName  = V3,
    Database    = V4,
    DomainStart = V5,
    DomainEnd   = V6
  )

cat("Consolidating domains...\n")
domains_consolidated <- cluster_domains(domains)

cat("Writing to", output_file, "\n")
fwrite(domains_consolidated, output_file)

cat("Done!\n")