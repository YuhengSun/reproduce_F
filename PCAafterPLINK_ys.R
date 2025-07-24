library(tidyverse)
library(ggrepel)

pca <- read_table2("sparrows2.eigenvec", col_names = FALSE) # this is Francesco's eigenvectors, for comparison
pca2 <- read_table2("sparrows_ys.eigenvec", col_names = FALSE) # this is my eigenvectors
eigenval <- scan("sparrows2.eigenval") # Francesco's eigenvalues
eigenval2 <- scan("sparrows_ys.eigenval") # my eigenvalues
# So hereafter, everything without suffix is produced using Francesco's data
# everything with the suffix "2" is produced using my data

#metadata <- read_csv("/home/fq/Desktop/uni/metadata.csv", col_names = TRUE)
metadata <- read_csv("metadata.csv", col_names = TRUE)
metadata <- as.tibble(data.frame(metadata))

# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
pca2 <- pca2[,-1]
# set names
names(pca)[1] <- "ind"
names(pca2)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
names(pca2)[2:ncol(pca2)] <- paste0("PC", 1:(ncol(pca)-1))


# sort out the individual species and pops
# sex
# sex <- rep(NA, length(pca$ind))
# sex[grep("PDOM*M", pca$ind)] <- "male"
# sex[grep("PDOM*F", pca$ind)] <- "female"
# sex[grep("PDOM*U", pca$ind)] <- "unidentified"
 # location
loc <- rep(NA, length(pca$ind))
loc2 <- rep(NA, length(pca2$ind))
loc[grep("AUS", pca$ind)] <- "Australia"
loc[grep("AUS", pca$ind)] <- "Australia"
loc2[grep("NCL", pca2$ind)] <- "New Caledonia"
loc2[grep("AUS", pca2$ind)] <- "Australia"
# combine - if you want to plot each in different colours
#sex_loc <- paste0(sex, "_", loc)

#added by me to name the dots
# sample_name <- df$origin

pca <- as.tibble(data.frame(pca, loc))
pca2 <- as.tibble(data.frame(pca2, loc2))

# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
pve2 <- data.frame(PC = 1:20, pve = eigenval2/sum(eigenval2)*100)

#merge the 2 dataframes
df = merge(x=pca, y=metadata, by="ind")
df2 = merge(x=pca2, y=metadata, by="ind")

#add country code to "origin" column
df <- df %>%
  mutate(origin = paste0(origin, " (", country.code, ")"))
df2 <- df2 %>%
  mutate(origin = paste0(origin, " (", country.code, ")"))

# # Add a new column that indicates the continent
df <- df %>%
  mutate(region = case_when(
    country.code %in% c("AUS", "NCL") ~ "Oceania",
    country.code == "GBR" ~ "Europe",
    TRUE ~ "Europe"  # NA for all other cases
  ))
df2 <- df2 %>%
  mutate(region = case_when(
    country.code %in% c("AUS", "NCL") ~ "Oceania",
    country.code == "GBR" ~ "Europe",
    TRUE ~ "Europe"  # NA for all other cases
  ))

# for labeling the dots
# Calculate cluster centroids
 # centroids <- df %>%
 #   group_by(origin) %>%
 #   summarise(across(starts_with("PC"), mean))
 # centroids2 = merge(x=centroids, y=df, by="origin")
 # centroids2 <- centroids2 %>% mutate(hjust = c(-0.35, -0.35, -0.35, -0.35, -0.35, -0.35, -0.35, -0.35, -0.35, -0.35, -0.35, -0.35, -0.35, -0.35, -0.35))
 
# centroids <- df %>%
#   group_by(origin) %>%
#   summarise()

#plot the two merged dataframes (use this to plot metadata)
# plot pca
b <- ggplot(df, aes(PC1, PC2, col = origin, shape = region
                    )) + geom_point(size = 1.5) +
scale_colour_manual(values = c("red", "blue", "darkorange", "darkcyan", "gold","green4", "tan4", "indianred","skyblue","burlywood","darkseagreen","pink","yellow3","darkgoldenrod","tomato2")) +
theme_light() + #geom_text(data = centroids2, aes(PC1.x, PC2.x, label = origin), size = 3, hjust = -0.35) +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

b + guides(color = "none")

b

b2 <- ggplot(df2, aes(PC1, PC2, col = origin, shape = region
)) + geom_point(size = 1.5) +
  scale_colour_manual(values = c("red", "blue", "darkorange", "darkcyan", "gold","green4", "tan4", "indianred","skyblue","burlywood","darkseagreen","pink","yellow3","darkgoldenrod","tomato2")) +
  theme_light() + #geom_text(data = centroids2, aes(PC1.x, PC2.x, label = origin), size = 3, hjust = -0.35) +
  xlab(paste0("PC1 (", signif(pve2$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve2$pve[2], 3), "%)"))

b2 + guides(color = "none")

b2

#d <- ggplot(centroids2, aes(PC1.x, PC2.x, col = origin, label = origin)) + geom_point(size = 3)+
  geom_text_repel(data = centroids2, aes())
#d

c <- ggplot(df, aes(PC1, PC2, col = origin)) + geom_point(size = 1.5) + theme_light() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  ggtitle("title")

c

# make plot of variance explained by first 20 PC
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

# plot pca
b <- ggplot(pca, aes(PC1, PC2, col = loc, label = ind)) + geom_point(size = 2) 
b <- b + scale_colour_manual(values = c("red", "blue"))
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

