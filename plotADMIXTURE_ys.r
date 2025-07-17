# load packages
library(dplyr)
library(readr)

# read in the list file, keep only the first column (sample ID)
list_df <- read_table("sparrows_ys.list", col_names = FALSE)[1]

# read in metadata
meta_df <- read_csv("metadata.csv")%>%mutate(cluster = case_when(
    `country code` %in% c("NOR", "GBR", "NLD") ~ "Europe",
    origin %in% c("mount_isa", "richmond", "townsville") ~ "north_cluster",
    TRUE ~ "south_cluster"
  ))

# merge
info_df <- left_join(list_df, meta_df, by = c("X1" = "ind"))%>%
  select(X1, cluster)

# write out the info file
write.table(info_df, file = "sparrows_ys.info", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
