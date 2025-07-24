library(ggplot2)
library(dplyr)

# Read in missing rate data
miss <- read.table("out.imiss", header = TRUE)

# Read in depth data
depth <- read.table("depth.out.idepth", header = TRUE)

# Merge the two dataframes based on individual
merged <- merge(miss, depth, by = "INDV")%>%
  mutate(DEPTH_SCALED = MEAN_DEPTH / max(MEAN_DEPTH, na.rm = TRUE)) # Scale the depth for better visualisation

# Order and keep the order
# merged <- merged[order(merged$F_MISS, decreasing = TRUE), ]
# merged$INDV <- factor(merged$INDV, levels = merged$INDV)

# Plot, where bars are missing rates, and dots are depths
ggplot(merged, aes(x = INDV)) +
  geom_bar(aes(y = F_MISS), stat = "identity", fill = "#69b3a2") +
  geom_point(aes(y = DEPTH_SCALED), color = "red", size = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
  labs(x = "Individual", y = "Missing rate (bar) & scaled depth (red dot)",
       title = "Per-individual missing rate and sequencing depth")



