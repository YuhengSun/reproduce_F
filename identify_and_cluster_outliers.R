### identify putative outlier genes
#rm(list = ls())

#setwd("/home/fq/Desktop/uni/haplotype_statistics/cluster_outputs/")

library(tidyverse)

#filter by logpval

#myData <- read_csv("melbourne_noumea_xpehh_nochrZ.csv")
myData <- read_csv("../../xpEHH/OutlierF.csv")%>%filter(Chromosome!="chrZ")

#myData <- filter(myData, LOGPVALUE > 6)

myData <- arrange(myData, Chromosome, cumulative_position)

# cluster outlier regions - find max value within them
threshold <- 100000 #original script: 100000
chr_clust <- unlist(sapply(unique(myData$Chromosome), function(w){
  y <- filter(myData, Chromosome == w) %>% .$cumulative_position
  y <- cumsum(c(1, diff(y) > threshold))
  paste0(w, "_", y)
}))
names(chr_clust) <- NULL
myData$cluster <- factor(chr_clust)

myCluster <- myData %>% group_by(cluster) %>% summarise(start = min(cumulative_position), stop = max(cumulative_position), size = stop - start)

# # # this part, I added to fiddle around

myCluster <- arrange(myCluster, start)

# Adjust the start and stop positions
myCluster <- myCluster %>%
  mutate(start = start - 125000,
         stop = stop + 125000)

# # #

# identify highest peaks
myData <- myData %>% group_by(cluster) %>% filter(PVAL == max(PVAL))
#write.csv(myData, file = "adelaide_bhill_xpehh_nochrZ_outlier_clusters.csv", row.names = FALSE)
write.csv(myCluster, file = "melbourne_noumea_xpehh_nochrZ_clusters_+-125k.csv", row.names = FALSE)
