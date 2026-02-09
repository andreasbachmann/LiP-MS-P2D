library(data.table)
library(dplyr)

rbp2go <- fread("data/Search_Results_InterPro.csv")
custom_rnas <- fread("data/RNAbinding_IDs.tsv", header = FALSE, col.names = c("ID", "Name"))

# keep only RNA-binding domains
rbp2go <- rbp2go %>% filter(RNA_binding == TRUE)

# combine both
combined <- bind_rows(custom_rnas, rbp2go %>% select(ID, Name)) %>%
  distinct(ID, .keep_all = TRUE) %>%
  select(DomainName = Name)

# save output
fwrite(combined, "data/unique_rbd.tsv", sep = "\t")
