library(rehh)
library(tidyverse)

# Extract filename suffixes
suffix <- sub("^Adelaide_", "", list.files("Adelaide", pattern = "\\.vcf\\.gz$"))

### xpEHH

for (i in suffix) {
# read in data for each population
# Adelaide
Adelaide_hh <- data2haplohh(hap_file = paste0("Adelaide/Adelaide_", i),
                         polarize_vcf = FALSE)
# bactrianus
BrokenHill_hh <- data2haplohh(hap_file = paste0("BrokenHill/BrokenHill_", i),
                       polarize_vcf = FALSE)

# Filter on MAF
Adelaide_hh_f <- subset(Adelaide_hh, min_maf = 0.05)
BrokenHill_hh_f <- subset(BrokenHill_hh, min_maf = 0.05)

# Perform scans
Adelaide_scan <- scan_hh(Adelaide_hh_f, polarized = F)
BrokenHill_scan <- scan_hh(BrokenHill_hh_f, polarized = F)

# Perform iHS
Adelaide_ihs <- ihh2ihs(Adelaide_scan, freqbin = 1)
BrokenHill_ihs <- ihh2ihs(BrokenHill_scan, freqbin = 1)

# ggplot(Adelaide_ihs$ihs, aes(POSITION, IHS)) + geom_point()
# ggplot(Adelaide_ihs$ihs, aes(POSITION, LOGPVALUE)) + geom_point()

# Perform xp-ehh
Adelaide_BrokenHill <- ies2xpehh(Adelaide_scan, BrokenHill_scan,
                                 popname1 = "Adelaide", popname2 = "Broken Hill",
                                 include_freq = T)

save(Adelaide_BrokenHill, file = paste0("Adelaide_BrokenHill_", sub(".vcf.gz", "", i), ".RData"))
print(paste(sub(".vcf.gz", "", i), "saved."))
}

### Plot

# Arrange the suffix so the chromosomes will be in order
chr <- sub(".vcf.gz", "", suffix)
nums <- as.numeric(gsub("chr(\\d+).*", "\\1", chr))
chr <- c(
  chr[!is.na(nums)][order(nums[!is.na(nums)])],
  chr[is.na(nums)]
)
chr

# Initialize an empty data frame to store combined results
combined_data <- data.frame()
cumulative_position <- 0  # Initialize cumulative position
mydata <- 0 # Initialize the position where new chromosomes start

# Skip chr16 for some reason
for (i in chr[-16]) {
  # Load the data
  load(paste0("Adelaide_BrokenHill_", i, ".RData"))
  
  # Extract relevant data and add a cumulative position column
  xpehh_data <- select(Adelaide_BrokenHill, CHR, POSITION, `XPEHH_Adelaide_Broken Hill`, LOGPVALUE)
  xpehh_data$cumulative_position <- xpehh_data$POSITION + cumulative_position
  
  # Update cumulative position for the next chromosome
  cumulative_position <- max(xpehh_data$cumulative_position)
  
  # Combine data
  combined_data <- rbind(combined_data, xpehh_data)
  mydata <- c(mydata, cumulative_position)
  
  print(i)
}

summary(combined_data$LOGPVALUE)

tick_positions <- (mydata[-1] + mydata[-length(mydata)]) / 2
tick_positions

# create a color column based on the threshold
combined_data$outlier <-ifelse(combined_data$LOGPVALUE > 6, "outlier", "background")

#plot log p value of xpehh

plot1 <- ggplot(combined_data, aes(x = cumulative_position, y = `LOGPVALUE`, color = outlier)) +
  geom_point(size = 0.2) +
  scale_color_manual(values = c("outlier" = "dodgerblue", "background" = "lightcyan3")) +
  labs(title = "Adelaide - Broken Hill",
       x = " ",
       y = "log (p) xp-EHH" ) +
  theme_light() + #common_theme +a
  #scale_color_discrete(name = "Chromosome") 
  geom_vline(xintercept = mydata,
             linetype = "dashed", color = "grey1", linewidth = 0.1) +
  scale_x_continuous(breaks = tick_positions, labels = chr[-16]) +
  #  geom_hline(yintercept = 6, linewidth = 0.1, color = "black") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot1
