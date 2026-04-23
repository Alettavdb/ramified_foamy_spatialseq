# approach 7-8 

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(ggplot2)
  library(ggpubr)
})

args<-commandArgs(T)

if (length(args)<1) {
  stop("Usage: Rscript lym_neighbors.R [chip id]")
}

chipid<-args[1]
#knum <- args[2]

print(paste0("This is chip ", chipid))
allMeta<-read.csv("MS_data/202511_neighborhood_analysis/allMeta_lymType.csv", row.names = "X")

allMeta$lymType<-gsub("non-vascular", "parenchymal_lym", allMeta$lymType)
allMeta$lymType<-gsub("vascular", "perivascular_lym", allMeta$lymType)
allMeta$lymType<-gsub("perivascular_lym_cells", "vascular_cells", allMeta$lymType)


chipMeta<-allMeta[allMeta$id==chipid,]

table(chipMeta$lymType)


target_cellType<-"parenchymal_lym" # parenchymal_lym perivascular_lym

target_cells_df <- chipMeta %>% filter(lymType == target_cellType)
vascular_lym_df <- chipMeta %>% filter(lymType == "vascular_cells")
oligo_lym_df <- chipMeta %>% filter(lymType == "oligodendrocytes")
astro_lym_df <- chipMeta %>% filter(lymType == "astrocytes")
mg_lym_df <- chipMeta %>% filter(lymType == "microglia")

# --- T1: Perivascular Lymphocytes ---
data_t1_matrix <- as.matrix(vascular_lym_df[, c("coord_x", "coord_y")], rownames.force = FALSE)
mode(data_t1_matrix) <- "numeric"
data_t2_matrix <- as.matrix(oligo_lym_df[, c("coord_x", "coord_y")], rownames.force = FALSE)
mode(data_t2_matrix) <- "numeric"
data_t3_matrix <- as.matrix(astro_lym_df[, c("coord_x", "coord_y")], rownames.force = FALSE)
mode(data_t2_matrix) <- "numeric"
data_t4_matrix <- as.matrix(mg_lym_df[, c("coord_x", "coord_y")], rownames.force = FALSE)
mode(data_t2_matrix) <- "numeric"

# --- R: Vascular Cells (Query) ---
query_r_matrix <- as.matrix(target_cells_df[, c("coord_x", "coord_y")], rownames.force = FALSE)
mode(query_r_matrix) <- "numeric"

# --- Distance to Perivascular Lymphocytes (T1) ---
nabor_t1_result <- nabor::knn(
  data = data_t1_matrix,           # The reference pool
  query = query_r_matrix,          # The query points
  k = 1
)

# nabor stores the distances in a list element called 'nn.dist'
dist_t1 <- as.numeric(nabor_t1_result$nn.dist)

# --- Distance to Parenchymal Lymphocytes (T1) ---
nabor_t2_result <- nabor::knn(
  data = data_t2_matrix,           # The reference pool
  query = query_r_matrix,          # The query points
  k = 1
)
# nabor stores the distances in a list element called 'nn.dist'
dist_t2 <- as.numeric(nabor_t2_result$nn.dist)

nabor_t3_result <- nabor::knn(
  data = data_t3_matrix,           # The reference pool
  query = query_r_matrix,          # The query points
  k = 1
)
# nabor stores the distances in a list element called 'nn.dist'
dist_t3 <- as.numeric(nabor_t3_result$nn.dist)

nabor_t4_result <- nabor::knn(
  data = data_t4_matrix,           # The reference pool
  query = query_r_matrix,          # The query points
  k = 1
)
# nabor stores the distances in a list element called 'nn.dist'
dist_t4 <- as.numeric(nabor_t4_result$nn.dist)

# Combine the distances into a single data frame for comparison and visualization
distance_df <- data.frame(
  Lym_Cell_ID = rownames(target_cells_df),
  Distance_to_VascularCells = as.vector(dist_t1),
  Distance_to_Oligodendrocytes = as.vector(dist_t2),
  Distance_to_Astrocytes = as.vector(dist_t3),
  Distance_to_Microglia = as.vector(dist_t4)
)

dist_df2<-distance_df
dist_df2$location<-gsub("_lym","",target_cellType)
dist_df2$donor<-substring(dist_df2$Lym_Cell_ID, 1, 4)
write.table(dist_df2, paste0(chipid, "_distances.csv"), sep=",", quote=F, col.names=F, row.names=F)

#stop(paste0(chipid," Distances Exported."))

## plot the distances
# 1. Convert the data from wide to long format
distance_long <- distance_df %>%
  pivot_longer(
    cols = starts_with("Distance_to"), # Selects the two distance columns
    names_to = "Cell_Type",
    values_to = "Distance"
  ) %>%
  # Clean up the names for better plotting labels
  mutate(
    Cell_Type = gsub("Distance_to_", "", Cell_Type)
  )
head(as.data.frame(distance_long))

# Load the necessary library for post-hoc tests and visualization
library(ggpubr)
# Perform the Friedman Test
# The formula is: Response variable ~ Grouping variable | Blocking/Subject ID variable
friedman_result <- friedman.test(
  Distance ~ Cell_Type | Lym_Cell_ID,
  data = distance_long
)

print("--- Overall Friedman Test Result ---")
print(friedman_result)

# Dunn's test (post-hoc) with Bonferroni adjustment for multiple comparisons
dunn_result <- rstatix::dunn_test(
  Distance ~ Cell_Type,
  data = distance_long,
  p.adjust.method = "bonferroni"
)

print("--- Dunn's Post-Hoc Test Result (Pairwise) ---")
print(as.data.frame(dunn_result))


############## plot
pdf(paste0(chipid, "_lym2others_dist.pdf"), width=8, height=5)  ##################################################
my_comparison <- list(
  c("Perivascular", "Parenchymal")
)
ggplot(distance_long, aes(x = Cell_Type, y = Distance, fill = Cell_Type)) +
  # 1. Boxplot and Jitter
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), size = 0.5, alpha = 0.2) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white")+
# 2. Add the Overall Friedman Test P-value to the plot corner
  stat_compare_means(
    method = "kruskal.test", # Kruskal-Wallis is the unpaired equivalent, but ggpubr uses it for plotting multi-group summaries
    label.x.npc = "left",   # Position the p-value
    label.y.npc = "top",
    label = "p.format"
  ) +

# 4. Final Labels
  labs(
    title = paste0("Proximity of Perivascular Lym to Other Cell Types: ", chipid),
    subtitle = paste0("Overall Test: Friedman Test (P=", format.pval(friedman_result$p.value, digits=3), ")"),
    x = "Nearest Neighbor Type",
    y = "Distance",
    fill = "Cell Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none") # Hide redundant legend
dev.off()

#print(friedman_result)
