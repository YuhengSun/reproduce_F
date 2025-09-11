### identify putative outlier genes
#rm(list = ls())

#setwd("/home/fq/Desktop/uni/divstats_vcftools/fst/")

library(tidyverse)

#load fst values
#myData <- read_tsv("melbourne_noumea/fst_melbourne_noumea.windowed.weir.fst")
myData <- read_tsv("../../fst/fst_adelaide_brokenhill.windowed.weir.fst")
#remove scaffolds
myData <- myData[!grepl("^scaffold", myData$CHROM), ]
# Remove chrZ
myData <- myData[!grepl("chrZ", myData$CHROM), ]


#modify start and end columna so that it is sequential
prev_value <- 0

for (i in 1:nrow(myData)) {
  # Check if the text in the text column changes
  if (i > 1 && myData$CHROM[i] != myData$CHROM[i - 1]) {
    # Update the previous value
    prev_value <- myData$BIN_START[i - 1]
  }
  myData$BIN_START[i] <- myData$BIN_START[i] + prev_value
}

myData$BIN_END <- myData$BIN_START + 19999

# add bin mid column
myData$mid <- (myData$BIN_END + myData$BIN_START) / 2

# identify the 95% percentile
my_threshold1 <- quantile(myData$WEIGHTED_FST, 0.975, na.rm = T)
# make an outlier column in the data.frame
myData <- myData %>% mutate(outlier = ifelse(WEIGHTED_FST > my_threshold1, "outlier", "background"))

myData <- filter(myData, outlier == "outlier")
myData <- arrange(myData, CHROM, mid)

# cluster outlier regions - find max value within them
threshold <- 100000 #original script: 100000
chr_clust <- unlist(sapply(unique(myData$CHROM), function(w){
  y <- filter(myData, CHROM == w) %>% .$mid
  y <- cumsum(c(1, diff(y) > threshold))
  paste0(w, "_", y)
}))
names(chr_clust) <- NULL
myData$cluster <- factor(chr_clust)

myCluster <- myData %>% group_by(cluster) %>% summarise(start = min(BIN_START), stop = max(BIN_END), size = stop - start)

# # # this part, I added to fiddle around

myCluster <- arrange(myCluster, start)

# # Adjust the start and stop positions
# myCluster <- myCluster %>%
#   mutate(start = start - 125000,
#          stop = stop + 125000)

# # #

# identify highest peaks
myData <- myData %>% group_by(cluster) %>% filter(WEIGHTED_FST == max(WEIGHTED_FST))
write.csv(myCluster, file = "outlier_clusters/melbourne_noumea_top2,5%fst_nochrZ_clusters_20kb.csv", row.names = FALSE)
