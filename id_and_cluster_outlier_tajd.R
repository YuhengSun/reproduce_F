### identify putative outlier genes
rm(list = ls())

setwd("/home/fq/Desktop/uni/tajD")

library(tidyverse)

#load data

myData <- read_tsv("10kb/mtisa_10kb.Tajima.D")


# modify bin start column

prev_value <- 0

for (i in 1:nrow(myData)) {
  # Check if the text in the text column changes
  if (i > 1 && myData$CHROM[i] != myData$CHROM[i - 1]) { 
    #prev_value <- NA  # Reset previous value when text changes
    # Update the previous value
    prev_value <- myData$BIN_START[i - 1]
  }
  myData$BIN_START[i] <- myData$BIN_START[i] + prev_value
}

# add bin end column
myData$BIN_END <- myData$BIN_START + 9999
# add bin mid column
myData$mid <- (myData$BIN_END + myData$BIN_START) / 2

#remove scaffolds
taj.all <- myData[!grepl("^scaffold", myData$CHROM), ]
#remove chrZ
taj.all <- taj.all[!grepl("chrZ", taj.all$CHROM), ]

# identify the lowest 2.5%
my_threshold1 <- quantile(taj.all$TajimaD, 0.025, na.rm = T)
# make an outlier column in the data.frame
taj.all <- taj.all %>% mutate(outlier = ifelse(TajimaD < my_threshold1, "outlier", "background"))

# remove non-outliers
taj.all <- filter(taj.all, TajimaD < my_threshold1)

# cluster outlier regions - find max value within them
threshold <- 100000 #original script: 100000
chr_clust <- unlist(sapply(unique(taj.all$CHROM), function(w){
  y <- filter(taj.all, CHROM == w) %>% .$mid
  y <- cumsum(c(1, diff(y) > threshold))
  paste0(w, "_", y)
}))
names(chr_clust) <- NULL
taj.all$cluster <- factor(chr_clust)

myCluster <- taj.all %>% group_by(cluster) %>% summarise(start = min(BIN_START), stop = max(BIN_END), size = stop - start)

myCluster <- arrange(myCluster, start)

# # #

#write.csv(myData, file = "outlier_clusters/adelaide_bhill_xpehh_nochrZ_outlier_clusters.csv", row.names = FALSE)
write.csv(myCluster, file = "outlier_clusters_10kb/mtisa_tajimad_bottom2,5%_clusters.csv", row.names = FALSE)

### plot cumulative snp count
max_value <- max(myData$N_SNPS)
x_values <- 1:max_value
y_values <- sapply(x_values, function(x) sum(myData$N_SNPS >= x))
plot_data <- data.frame(x = x_values, y = y_values)
ggplot(plot_data, aes(x = x, y = y)) +
  geom_line() +
  labs(title = "Cumulative Count of N_snps Values",
       x = "N_snps Value",
       y = "Number of Rows with N_snps >= x") +
  theme_minimal()