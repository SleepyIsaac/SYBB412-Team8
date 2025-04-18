library(tibble)
library(Seurat)
library(tidyverse)

# Load individual Seurat objects
TN_seurat <- LoadSeuratRds("../data/TN.rds")
PD_seurat <- LoadSeuratRds("../data/PD.rds")
RD_seurat <- LoadSeuratRds("../data/RD.rds")

# Merge datasets and save merged object
merged <- merge(TN_seurat, y = list(RD_seurat, PD_seurat), add.cell.ids = c("TN", "RD", "PD"))
saveRDS(merged, "../data/merged.rds")

# Extract Smart-seq2 count matrices for each group
counts_TN <- merged@assays$RNA@layers[["counts.SmartSeq2_TN"]]
counts_RD <- merged@assays$RNA@layers[["counts.SmartSeq2_RD"]]
counts_PD <- merged@assays$RNA@layers[["counts.SmartSeq2_PD"]]

# Rename columns to match cell barcodes
colnames_TN <- grep("^TN_", colnames(merged), value = TRUE)
colnames_RD <- grep("^RD_", colnames(merged), value = TRUE)
colnames_PD <- grep("^PD_", colnames(merged), value = TRUE)

colnames(counts_TN) <- colnames_TN
colnames(counts_RD) <- colnames_RD
colnames(counts_PD) <- colnames_PD

# Combine count matrices
expr_mat <- cbind(counts_TN, counts_RD, counts_PD)
rownames(expr_mat) <- rownames(merged@assays$RNA)

# Convert to data frame and export count matrix
counts_df <- as.data.frame(as.matrix(expr_mat))
counts_df <- rownames_to_column(counts_df, var = "gene")
write.table(counts_df, file = "../data/counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Prepare and export cell metadata (for CellPhoneDB or Monocle)
meta_cpdb <- data.frame(
  cell = colnames(expr_mat),
  cell_type = merged$cell_type[colnames(expr_mat)]
)
write.table(meta_cpdb, file = "../data/meta.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
