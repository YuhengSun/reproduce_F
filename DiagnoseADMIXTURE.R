### diagnose ADMIXTURE result

library(tidyverse)

# Read in Q files for K=2 and K=3
q2 <- read.table("sparrows_ys.2.Q", header = F)
q3 <- read.table("sparrows_ys.3.Q", header = F)
# Individual ID:
ind <- read.table("sparrows_ys.fam")$V1
# Cluster info
cluster <- read.table("sparrows_ys.info", col.names = c("ind", "cluster"))
# Integrate the dataframes
q2_full <- cbind(ind, q2)%>%left_join(cluster)

df <- q2_full

# Set cluster as a factor with order, and arrange it to make it easier to compare with Francesco's result
df$cluster <- factor(df$cluster, levels = c("north_cluster", "south_cluster", "Europe"))

# Arrange by cluster
df <- df %>% arrange(cluster)

# Assign each individual an order, and this will be used as the x axis
df$ind_order <- factor(df$ind, levels = df$ind)

# Convert into long format for ggplot to plot stacked bars
df_long <- df %>%
  pivot_longer(cols = c(V1, V2), names_to = "Ancestry", values_to = "Proportion")%>%
  mutate(Ancestry = factor(Ancestry, levels = c("V2", "V1")))

# Plot
ggplot(df_long, aes(x = ind_order, y = Proportion, fill = Ancestry)) +
  geom_bar(stat = "identity", width = 1) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank()
  ) +
  facet_grid(~ cluster, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("#2c7bb6", "#d7191c")) + 
  xlab("Individuals") +
  ylab("Ancestry Proportion") +
  ggtitle("ADMIXTURE Barplot (K = 2)")

# Integrate the dataframes for K=3
q3_full <- cbind(ind, q3)%>%left_join(cluster)

# Do the same
df <- q3_full
df$cluster <- factor(df$cluster, levels = c("north_cluster", "south_cluster", "Europe"))
df <- df %>% arrange(cluster)
df$ind_order <- factor(df$ind, levels = df$ind)
df_long <- df %>%
  pivot_longer(cols = c(V1, V2, V3), names_to = "Ancestry", values_to = "Proportion") %>%
  mutate(Ancestry = factor(Ancestry, levels = c("V1", "V2", "V3")))
ggplot(df_long, aes(x = ind_order, y = Proportion, fill = Ancestry)) +
  geom_bar(stat = "identity", width = 1) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank()
  ) +
  facet_grid(~ cluster, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("#4daf4a", "#377eb8", "#e41a1c")) + 
  xlab("Individuals") +
  ylab("Ancestry Proportion") +
  ggtitle("ADMIXTURE Barplot (K = 3)")
