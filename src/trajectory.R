library(monocle3)
library(Matrix)
library(ggplot2)

# Load count matrix and metadata
counts <- read.delim("../data/counts.tsv", row.names = 1, check.names = FALSE)
meta <- read.delim("../data/meta.tsv")

# Convert expression matrix to sparse matrix format
expr_matrix <- as(as.matrix(counts), "dgCMatrix")

# Prepare cell and gene metadata
cell_metadata <- meta
rownames(cell_metadata) <- cell_metadata$cell
gene_metadata <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))

# Create a monocle3 cell_data_set (CDS) object
cds <- new_cell_data_set(expr_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata)

# Preprocess and reduce dimensions
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds)

# Cluster cells and learn trajectory graph
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

# Order cells along pseudotime based on UMAP
cds <- order_cells(cds, reduction_method = "UMAP")

# Select genes of interest for visualization
selected_genes <- c("TP53", "EGFR", "PDCD1", "CD274",
                    "CEACAM5", "IFNG", "MMP2", "SOX2", "ALDH1A1")

# Subset CDS to selected genes
cds_subset <- cds[rowData(cds)$gene_short_name %in% selected_genes, ]

# Plot gene expression over pseudotime
pseudotime_expression_plot <- plot_genes_in_pseudotime(cds_subset,
                                                       min_expr = 0.1,
                                                       color_cells_by = "pseudotime")

# Save plot
ggsave("../figures/pseudotime_expression.png", pseudotime_expression_plot, dpi = 300)
