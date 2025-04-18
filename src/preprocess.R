library(Seurat)
library(dplyr)
library(biomaRt)
library(ggplot2)
library(grid)
library(patchwork)
library(cowplot)

set.seed(412)

# Function to load, filter and create Seurat object for one group (e.g. TN, PD, RD)
create_filtered_seurat <- function(group_label) {
  file_name <- paste0("../data/", group_label, "_counts.txt")
  
  # Load count matrix
  counts <- read.delim(file_name, comment.char = "#", row.names = 1)
  counts[is.na(counts)] <- 0
  gene_id <- rownames(counts)
  counts <- apply(counts, 2, as.numeric)
  rownames(counts) <- gene_id
  counts <- counts[, -1]  # Remove first column if needed
  counts <- as.matrix(counts)
  
  # Clean column names
  colnames(counts)[5:ncol(counts)] <- sub(".*(SRR[0-9]+).*", "\\1", colnames(counts)[5:ncol(counts)])
  
  # Load metadata
  metadata <- read.csv("../data/metadata.csv")
  metadata <- metadata[metadata$treatement.timepoint == group_label, ]
  sample_ids <- metadata$sample_id
  rownames(metadata) <- metadata$run_id
  
  # Convert Ensembl ID to HGNC symbol
  ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
  gene_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                    filters = "ensembl_gene_id",
                    values = rownames(counts),
                    mart = ensembl)
  gene_map <- gene_map[gene_map$hgnc_symbol != "", ]
  gene_map <- gene_map[!duplicated(gene_map$ensembl_gene_id), ]
  
  # Keep only matched genes and set gene symbols
  matched_ids <- gene_map$ensembl_gene_id
  counts <- counts[rownames(counts) %in% matched_ids, ]
  gene_symbols <- gene_map$hgnc_symbol[match(rownames(counts), gene_map$ensembl_gene_id)]
  rownames(counts) <- gene_symbols
  counts <- counts[!duplicated(rownames(counts)), ]
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = counts, project = paste0("SmartSeq2_", group_label))
  
  # Get mitochondrial and ribosomal gene names
  mt_genes <- getBM(attributes = c("hgnc_symbol"), filters = "chromosome_name", values = "MT", mart = ensembl)$hgnc_symbol
  #ribo_genes <- getBM(attributes = c("hgnc_symbol"), filters = "chromosome_name", values = c(paste0("RPL", 1:50), paste0("RPS", 1:50)), mart = ensembl)$hgnc_symbol
  mt_genes <- intersect(mt_genes, rownames(seurat_obj))
  #ribo_genes <- intersect(ribo_genes, rownames(seurat_obj))
  
  ribo_symbols <- c(paste0("RPL", 1:50), paste0("RPS", 1:50))
  ribo_genes <- getBM(attributes = c("hgnc_symbol"),
                      filters = "hgnc_symbol",
                      values = ribo_symbols,
                      mart = ensembl)$hgnc_symbol
  
  
  # Calculate QC metrics
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, features = mt_genes)
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, features = ribo_genes)
  
  # Filter cells based on QC
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 500 & nFeature_RNA < 12000 & percent.mt < 15 & percent.ribo < 10)
  
  # Add metadata
  seurat_obj <- AddMetaData(seurat_obj, metadata = metadata)
  
  vln_plot <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)
  
  ggsave(paste0("../figures/CellQC_", group_label, ".png"), vln_plot, width = 16, height = 6, dpi = 300)
  
  return(list(seurat_obj = seurat_obj, group_label = group_label))
}

# Function to run clustering
find_clusters <- function(seurat, group_label, genes) {
  seurat <- NormalizeData(seurat)
  seurat <- FindVariableFeatures(seurat)
  seurat <- ScaleData(seurat)
  seurat <- RunPCA(seurat)
  seurat <- RunUMAP(seurat, dims = 1:30)
  #seurat_umap <- seurat
  seurat <- FindNeighbors(seurat)
  seurat <- FindClusters(seurat, resolution = 0.5)
  UMAP_plot <- DimPlot(seurat, reduction = "umap", label = TRUE, group.by = "ident")
  ggsave(paste0("../figures/UMAP_", group_label, ".png"), UMAP_plot, width = 8, height = 6, dpi = 300)
  
  feature_plot <- FeaturePlot(seurat, features = genes)
  ggsave(paste0("../figures/fea_", group_label, ".png"), feature_plot, width = 8, height = 6, dpi = 300)
  
  return(seurat)
}

# Load, filter, cluster for each group

genes = c("TP53", "EGFR", "PDCD1", "CD274",
          "CEACAM5", "IFNG", "MMP2", "SOX2", "ALDH1A1")

TN_seurat <- create_filtered_seurat("TN")
TN_seurat <- find_clusters(TN_seurat$seurat_obj, TN_seurat$group_label, genes)
PD_seurat <- create_filtered_seurat("PD")
PD_seurat <- find_clusters(PD_seurat$seurat_obj, PD_seurat$group_label, genes)
RD_seurat <- create_filtered_seurat("RD")
RD_seurat <- find_clusters(RD_seurat$seurat_obj, RD_seurat$group_label, genes)

# Define cluster annotations
cluster_annot <- c(
  "0" = "T cell", 
  "1" = "B cell", 
  "2" = "Macrophage", 
  "3" = "cDC",
  "4" = "Mast cell", 
  "5" = "Fibroblast", 
  "6" = "Endothelial", 
  "7" = "Epithelial",
  "8" = "pDC", 
  "9" = "Melanocyte"
)



# Add cell type annotation to each Seurat object
cluster_labels <- cluster_annot[as.character(TN_seurat$seurat_clusters)]
names(cluster_labels) <- colnames(TN_seurat)
TN_seurat <- AddMetaData(TN_seurat, metadata = cluster_labels, col.name = "cell_type")

cluster_labels <- cluster_annot[as.character(PD_seurat$seurat_clusters)]
names(cluster_labels) <- colnames(PD_seurat)
PD_seurat <- AddMetaData(PD_seurat, metadata = cluster_labels, col.name = "cell_type")

cluster_labels <- cluster_annot[as.character(RD_seurat$seurat_clusters)]
names(cluster_labels) <- colnames(RD_seurat)
RD_seurat <- AddMetaData(RD_seurat, metadata = cluster_labels, col.name = "cell_type")


# Save processed Seurat objects
saveRDS(TN_seurat, file = "../data/TN.rds")
saveRDS(PD_seurat, file = "../data/PD.rds")
saveRDS(RD_seurat, file = "../data/RD.rds")

# Define marker genes for dot plot and heatmap
marker_list <- list(
  "Tcell" = c("CD2", "CD3D", "CD3E", "CD3G"),
  "Bcell" = c("IGLL5", "MS4A1", "CD79A", "PAX5", "SDC1"),
  "Macrophage" = c("MARCO", "CD68", "CSF1R", "FOLR2", "C1QA", "APOE"),
  "Dendritic" = c("CD1C", "FCER1A", "NDRG2"),
  "Mast" = c("TPSAB1", "TPSB2", "CMA1"),
  "Neutrophils" = c("CXCR1", "CXCR2", "FCGR3B"),
  "Fibroblasts" = c("FAP", "THY1", "ACTA2"),
  "Endothelial" = c("PECAM1", "CD34", "VWF"),
  "Epithelial" = c("EPCAM", "KRT19", "SFN"),
  "pDCs" = c("CLEC4C", "LILRA4", "DNASE1L3"),
  "Melanocytes" = c("PMEL", "MLANA"),
  "Housekeeping" = c("ACTB", "GAPDH", "MALAT1")
)

features <- unique(unlist(marker_list))

# Dot plots for each condition
p_TN_dot <- DotPlot(TN_seurat, features = features, group.by = "seurat_clusters") +
  RotatedAxis() + ggtitle("TN N=663") +
  theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8))

p_PD_dot <- DotPlot(PD_seurat, features = features, group.by = "seurat_clusters") +
  RotatedAxis() + ggtitle("PD N=702") +
  theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), legend.position = "none")

p_RD_dot <- DotPlot(RD_seurat, features = features, group.by = "seurat_clusters") +
  RotatedAxis() + ggtitle("RD N=707") +
  theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), legend.position = "none")

# Add axis labels
x_label <- textGrob("Features", gp = gpar(fontsize = 14))
y_label <- textGrob("Cluster Identity", rot = 90, gp = gpar(fontsize = 14))

# Combine dot plots
dot_plot <- (p_TN_dot / p_PD_dot / p_RD_dot) + plot_layout(guides = "collect") &
  theme(axis.title = element_blank(), legend.position = "right")

dot_plot <- (
  plot_spacer() + plot_spacer() + plot_spacer() +
    wrap_elements(y_label) + dot_plot + wrap_elements(plot_spacer()) +
    plot_spacer() + wrap_elements(x_label) + plot_spacer()
) +
  plot_layout(ncol = 3, widths = c(0.05, 0.9, 0.05), heights = c(0.05, 0.9, 0.05))

ggsave("../figures/dot.png", dot_plot, width = 16, height = 12, dpi = 300)

# Heatmap colors
cluster_colors <- c(
  "0" = "#E41A1C", "1" = "#377EB8", "2" = "#4DAF4A", "3" = "#984EA3",
  "4" = "#FF7F00", "5" = "#FFFF33", "6" = "#A65628", "7" = "#F781BF",
  "8" = "#999999", "9" = "#66C2A5"
)

annot_colors <- setNames(cluster_colors[names(cluster_annot)], cluster_annot)

# Heatmap for each group
TN_heat <- DoHeatmap(TN_seurat, features = features, group.colors = cluster_colors) +
  ggtitle("TN N=663") +
  theme(strip.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 8), 
        legend.position = "none")

PD_heat <- DoHeatmap(PD_seurat, features = features, group.colors = cluster_colors) +
  ggtitle("PD N=702") +
  theme(axis.text.y = element_text(size = 8), 
        legend.position = "none")

RD_heat <- DoHeatmap(RD_seurat, features = features, group.colors = cluster_colors) +
  ggtitle("RD N=707") +
  theme(axis.text.y = element_text(size = 8), 
        legend.position = "right")

# Create custom legend
legend_df <- data.frame(
  cluster = factor(names(cluster_annot), levels = names(cluster_annot)),
  cell_type = cluster_annot
)


legend_plot <- ggplot(legend_df, aes(x = 1, fill = cell_type)) +
  geom_bar() +
  scale_fill_manual(
    values = annot_colors,
    name = "Cell type"
  ) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

# Combine heatmaps
heatmap_plot <- (TN_heat / PD_heat / RD_heat) + plot_layout(guides = "collect") &
  theme(legend.position = "right")

# Extract legend only
legend <- get_legend(legend_plot)

# Add legend to the right of heatmap
heatmap_plot <- plot_grid(
  heatmap_plot, legend,
  ncol = 2,
  rel_widths = c(4, 1)
)

ggsave("../figures/heatmap.png", heatmap_plot, width = 16, height = 12, dpi = 300)
