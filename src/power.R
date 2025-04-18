library(RNASeqPower)
library(dplyr)
library(tidyr)
library(ggplot2)

deg_results <- readRDS("../data/deg_results.rds")
effect_df <- deg_results %>%
  filter(!is.na(logFC)) %>%
  group_by(cluster_id, contrast) %>%
  summarise(effect = mean(abs(logFC)), .groups = "drop")


# Load individual Seurat objects
TN_seurat <- LoadSeuratRds("../data/TN.rds")
PD_seurat <- LoadSeuratRds("../data/PD.rds")
RD_seurat <- LoadSeuratRds("../data/RD.rds")

# Merge Seurat objects into one
merged <- merge(
  TN_seurat,
  y = list(RD_seurat, PD_seurat),
  add.cell.ids = c("TN", "RD", "PD")
)

# Set group/sample/cell type information in metadata
merged$group <- merged@meta.data$treatement.timepoint
merged$sample_id <- merged@meta.data$run_id
merged$cell_type <- merged@meta.data$cell_type

# Prepare gene names and extract counts matrix from RNA assay
gene_names <- rownames(merged@assays$RNA)
counts_TN <- merged@assays$RNA@layers[["counts.SmartSeq2_TN"]]
counts_RD <- merged@assays$RNA@layers[["counts.SmartSeq2_RD"]]
counts_PD <- merged@assays$RNA@layers[["counts.SmartSeq2_PD"]]
expr_mat <- cbind(counts_TN, counts_RD, counts_PD)
rownames(expr_mat) <- gene_names
meta <- merged@meta.data

# Create SingleCellExperiment (SCE) object
sce <- SingleCellExperiment(
  assays = list(counts = expr_mat),
  colData = meta
)

# Set required column names for muscat
colData(sce)$cluster_id <- sce$cell_type
colData(sce)$group_id <- as.character(sce$group)
colData(sce)$sample_id <- sce$sample_id

# Format SCE object for muscat
sce <- prepSCE(
  sce,
  kid = "cluster_id",  # cell type
  gid = "group_id",    # condition/group
  sid = "sample_id"    # sample ID
)

# Generate pseudobulk data
pb <- aggregateData(sce)

meta <- as.data.frame(colData(pb))

effect_df <- effect_df %>%
  mutate(group1 = sub(" - .*", "", contrast),
         group2 = sub(".* - ", "", contrast))

get_n_per_group <- function(cluster, g1, g2, meta) {
  sub_meta <- meta %>% filter(cell_type == cluster & group_id %in% c(g1, g2))
  n1 <- n_distinct(sub_meta$run_id[sub_meta$group_id == g1])
  n2 <- n_distinct(sub_meta$run_id[sub_meta$group_id == g2])
  return(min(n1, n2))
}

effect_df$n <- mapply(get_n_per_group,
                      cluster = effect_df$cluster_id,
                      g1 = effect_df$group1,
                      g2 = effect_df$group2,
                      MoreArgs = list(meta = meta))

estimate_power_batch <- function(effect_df, depth_million = 10, cv = 0.4, alpha = 0.05) {

  effect_df$power <- NA_real_
  
  for (i in 1:nrow(effect_df)) {
    row <- effect_df[i, ]
    
    if (row$n < 2 || is.na(row$effect)) {
      next
    }
    
    power_val <- tryCatch({
      rnapower(
        depth = depth_million,
        n = row$n,
        cv = cv,
        effect = row$effect,
        alpha = alpha
      )
    }, error = function(e) NA_real_)
    
    effect_df$power[i] <- power_val
  }
  
  return(effect_df)
}

power_df <- estimate_power_batch(effect_df, depth_million = 10, cv = 0.4)



power_df$label <- paste(power_df$cluster_id, power_df$contrast, sep = " / ")

power_df$label <- factor(power_df$label, levels = power_df$label[order(power_df$power)])

power_plot <- ggplot(power_df, aes(x = label, y = power, fill = power >= 0.8)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("red", "steelblue"), labels = c("< 0.8", "≥ 0.8")) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray40") +
  labs(
    title = "Estimated Power per Cluster / Contrast",
    x = "Cluster / Contrast",
    y = "Power",
    fill = "Power"
  ) +
  theme_minimal(base_size = 13)

ggsave("../figures/power_plot.png", power_plot, dpi=300)

high_power_pairs <- power_df %>%
  filter(power >= 0.8) %>%
  dplyr::select(cluster_id, contrast)

deg_results_filtered <- deg_results %>%
  inner_join(high_power_pairs, by = c("cluster_id", "contrast"))

saveRDS(deg_results_filtered, file = "../data/deg_results_high_power.rds")
