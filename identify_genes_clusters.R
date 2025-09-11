library(tidyverse)
library(stringr)
library(gtools)
library(dplyr)
library(tidyr)

rm(list = ls())

setwd("/home/fq/Desktop/uni")

#get gff file and name the columns
gff <- read_tsv("/home/fq/work/house_sparrow.gff", , col_names = FALSE)
colnames(gff) <- c("chr", "source", "feature", "start", "end", "score",
                   "strand", "frame", "attribute")

# select genes only
new_gff <- gff %>% filter(feature == "gene")
#remove scaffolds
new_gff <- new_gff %>% filter(!str_detect(chr, "scaffold"))
# Remove chrZ
new_gff <- new_gff %>% filter(!str_detect(chr, "chrZ"))
# make a gene mid point variable (temporary workaround - to be added in the for loop)
new_gff <- new_gff %>% mutate(mid = start + (end-start)/2)

#order all the chromosomes based on the value in the "mid" field
new_gff <- new_gff %>% group_by(chr) %>%   arrange(mid, .by_group = TRUE)

#order so that chr9 is before chr10
new_gff <- new_gff[mixedorder(new_gff$chr), ]
#order so that chr1a is first
# Filter rows where chr is "chr1A" and those that are not
chr1A_rows <- new_gff %>% filter(chr == "chr1A")
other_rows <- new_gff %>% filter(chr != "chr1A")
# Combine the filtered rows with chr1A_rows at the top
new_gff <- bind_rows(chr1A_rows, other_rows)

# make a gene mid point variable
prev_value <- 0
for (i in 1:nrow(new_gff)) {
  # Check if the text in the text column changes
  if (i > 1 && new_gff$chr[i] != new_gff$chr[i - 1]) { 
    # Update the previous value
    prev_value <- new_gff$mid[i - 1]
  }
  new_gff$mid[i] <- new_gff$mid[i] + prev_value
}

# load fst/xpehh/tajd data
w2000 <- read.csv("haplotype_statistics/cluster_outputs/adelaide_bhill_xpehh_nochrZ_clusters_+-125k.csv")

# add to w2000 the contents of new_gff$attribute if mid falls between start and stopp
# also add a count of the genes falling in each cluster

w2000 <- w2000 %>%
  rowwise() %>%  # Process each row independently
  mutate(
    # Concatenated attributes
    genes = ifelse(
      any(new_gff$mid >= start & new_gff$mid <= stop),
      paste(new_gff$attribute[new_gff$mid >= start & new_gff$mid <= stop], collapse = "; "),
      NA_character_
    ),
    # Count of 'mid' values falling within the range
    mid_count = sum(new_gff$mid >= start & new_gff$mid <= stop)
  ) %>%
  ungroup()
  
  # Write the output to a text file
write.csv(w2000, file = "haplotype_statistics/genes_in_clusters_outputs/Xcandidate_genes_clusters_xpehh_adelaide_bhill_nochrZ.csv", row.names = FALSE)

# extract only the genes in a list, for go analysis
gene_list <- str_extract(w2000$genes, "(?<=ID=)[^;]+")
gene_list <- na.omit(gene_list)

# write to text file
writeLines(gene_list, "haplotype_statistics/genes_in_clusters_outputs/XonlyID_xpehh_nochrZ_adelaide_bhill.list")

