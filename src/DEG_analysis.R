library(purrr)
library(dplyr)
library(Seurat)
library(muscat)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(tidyverse)
library(edgeR)

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

# Function: run edgeR on pseudobulk data for each cluster
run_edgeR_pb_clusterwise <- function(pb, contrasts_list) {
  meta <- as.data.frame(colData(pb))
  clusters <- names(assays(pb))
  results <- list()
  
  for (clust in clusters) {
    message("Processing cluster: ", clust)
    
    counts <- assay(pb, clust)
    if (is.null(rownames(counts))) stop("Counts matrix has no gene names.")
    
    # Filter out samples with no counts
    lib_sizes <- colSums(counts)
    counts <- counts[, lib_sizes > 0, drop = FALSE]
    meta_sub <- meta[colnames(counts), , drop = FALSE]
    
    # Check if all required groups are present
    valid_groups <- names(which(table(meta_sub$group_id) >= 2))
    if (length(valid_groups) < 2 || !all(c("TN", "PD", "RD") %in% valid_groups)) {
      message("Skipping cluster ", clust, ": missing TN or other groups.")
      next
    }
    
    # Build design matrix and run edgeR pipeline
    dge <- DGEList(counts = counts)
    dge$samples$group <- meta_sub$group_id
    meta_sub$group_id <- factor(meta_sub$group_id, levels = c("PD", "RD", "TN"))
    design <- model.matrix(~ 0 + group_id, data = meta_sub)
    colnames(design) <- gsub("group_id", "", colnames(design))
    dge <- calcNormFactors(dge)
    dge <- estimateDisp(dge, design)
    fit <- glmQLFit(dge, design)
    
    # Run edgeR for each contrast
    for (cname in names(contrasts_list)) {
      g1 <- contrasts_list[[cname]][1]
      g2 <- contrasts_list[[cname]][2]
      
      if (!(g1 %in% meta_sub$group_id) || !(g2 %in% meta_sub$group_id)) {
        message("Skipping contrast ", cname, " in cluster ", clust)
        next
      }
      
      res <- tryCatch({
        contrast_vec <- makeContrasts(
          contrasts = paste0(g1, " - ", g2),
          levels = design
        )
        test <- glmQLFTest(fit, contrast = contrast_vec)
        top <- topTags(test, n = Inf, sort.by = "none")$table
        
        top$gene <- rownames(counts)[match(rownames(top), rownames(counts))]
        top$cluster_id <- clust
        top$contrast <- paste(g1, "-", g2)
        top$logCPM <- rowMeans(edgeR::cpm(dge, log = TRUE)[rownames(top), ])
        top
      }, error = function(e) {
        message("Error in cluster ", clust, " contrast ", cname)
        NULL
      })
      
      if (!is.null(res)) {
        results[[paste(clust, cname, sep = "_")]] <- res
      }
    }
  }
  
  # Combine all results into one data frame
  bind_rows(results)
}

# Define contrasts to test
my_contrasts <- list(
  "PD_RD" = c("PD", "RD"),
  "TN_PD" = c("TN", "PD"),
  "TN_RD" = c("TN", "RD")
)

# Run DEG analysis
deg_results <- run_edgeR_pb_clusterwise(pb, my_contrasts)

# Save all DEG results
saveRDS(deg_results, file = "../data/deg_results.rds")

# Create a combined ranked list for GSEA
gsea_all <- deg_results %>%
  filter(!is.na(PValue), !is.na(logFC)) %>%
  mutate(rank_metric = logFC * -log10(PValue)) %>%
  arrange(desc(rank_metric)) %>%
  dplyr::select(gene, rank_metric)

# Write combined GSEA .rnk file
write.table(gsea_all, file = "../GSEA_rnk_files/GSEA_ranked_list.rnk",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Write per-cluster and contrast GSEA .rnk files
# deg_results %>%
#   filter(!is.na(PValue), !is.na(logFC)) %>%
#   mutate(rank_metric = logFC * -log10(PValue)) %>%
#   group_by(cluster_id, contrast) %>%
#   group_walk(~ {
#     contrast_clean <- gsub(" - ", "_", .y$contrast)
#     file_name <- paste0("../GSEA_rnk_files/cluster", "_", .y$cluster_id, "_", .y$contrast, ".rnk")
#     rnk_tbl <- .x %>% dplyr::select(gene, rank_metric) %>% arrange(desc(rank_metric))
#     write.table(rnk_tbl, file = file_name, sep = "\t",
#                 row.names = FALSE, col.names = FALSE, quote = FALSE)
#   })

# Add significance and -log10 FDR for volcano plot
deg_results$Significant <- with(deg_results, FDR < 0.05 & abs(logFC) > 1)
deg_results$neglogFDR <- -log10(deg_results$FDR)

# Create volcano plot faceted by cluster and contrast
volcano_plot <- ggplot(deg_results, aes(x = logFC, y = neglogFDR, color = Significant)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~ paste(cluster_id, contrast), scales = "free") +
  theme_minimal() +
  scale_color_manual(values = c("grey", "red")) +
  labs(y = "-log10 FDR")

# Save volcano plot image
ggsave("../figures/volcano_deg.png", volcano_plot, width = 16, height = 12, dpi = 300)

# Extract top 10 DEGs per group
top_deg <- deg_results %>%
  filter(FDR < 0.05) %>%
  group_by(cluster_id, contrast) %>%
  arrange(FDR) %>%
  slice_head(n = 10)

# Save DEG results
write.csv(deg_results, file = "../DEG_results/Full_DEG_Results.csv", row.names = FALSE)
write.csv(top_deg, file = "../DEG_results/Top10_DEGs.csv", row.names = FALSE)
write.csv(gsea_all, "../DEG_results/GSEA_ranked_combined.csv", row.names = FALSE)

# Save ranked DEG tables per cluster and contrast
deg_results %>%
  filter(!is.na(PValue), !is.na(logFC)) %>%
  mutate(rank_metric = logFC * -log10(PValue)) %>%
  group_by(cluster_id, contrast) %>%
  group_walk(~ {
    contrast_clean <- gsub(" - ", "_", .y$contrast)
    file_name <- paste0("../DEG_results/cluster", "_", .y$cluster_id, "_", contrast_clean, ".csv")
    ranked_tbl <- .x %>% dplyr::select(gene, rank_metric) %>% arrange(desc(rank_metric))
    write.csv(ranked_tbl, file = file_name, row.names = FALSE)
  })
