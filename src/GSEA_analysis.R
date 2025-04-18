library(fgsea)
library(msigdbr)
library(dplyr)
library(ggplot2)

# Load MSigDB gene sets for Biological Process (BP) category
pathways <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP") %>%
  split(x = .$gene_symbol, f = .$gs_name)

# Read in differential expression results
deg_results <- readRDS("../data/deg_results_high_power.rds")

# Define a function to run fgsea in batch mode for multiple contrasts and clusters
run_batch_fgsea <- function(deg_df, pathways, top_n = 30) {
  gsea_results <- list()  # Create an empty list to store GSEA results
  
  # Preprocess DEG data: remove NA values and calculate ranking metric
  deg_df <- deg_df %>%
    filter(!is.na(PValue), !is.na(logFC)) %>%
    mutate(rank_metric = logFC * -log10(PValue)) %>%
    group_by(cluster_id, contrast)
  
  # Split data by group (cluster + contrast)
  groups <- group_split(deg_df, .keep = TRUE)
  keys <- group_keys(deg_df)
  
  # Loop through each group
  for (i in seq_along(groups)) {
    .x <- groups[[i]]
    clust <- keys$cluster_id[[i]]
    cont <- keys$contrast[[i]]
    
    # Create ranked gene list
    rnk_tbl <- .x %>% dplyr::select(gene, rank_metric) %>% arrange(desc(rank_metric))
    rnk_vector <- deframe(rnk_tbl)
    
    # Check how many genes overlap with pathway gene sets
    matched_genes <- intersect(names(rnk_vector), unlist(pathways))
    current <- paste0(clust, "_", gsub(" - ", "_", cont))
    
    # Skip if too few genes matched
    if (length(matched_genes) < 5) {
      message("Skipping ", current, ": too few matched genes")
      next
    }
    
    # Run fgsea and handle potential errors
    fgsea_res <- tryCatch({
      fgsea(pathways = pathways, stats = rnk_vector, eps = 0.0)
    }, error = function(e) {
      message("Error in ", current, ": ", e$message)
      return(NULL)
    })
    
    # Skip if no results returned
    if (is.null(fgsea_res) || nrow(fgsea_res) == 0) {
      message("Skipping ", current, ": no results.")
      next
    }
    
    # Store result
    gsea_results[[current]] <- fgsea_res
    
    # Save result as CSV
    fgsea_res <- fgsea_res %>% arrange(padj)
    fgsea_res$leadingEdge <- sapply(fgsea_res$leadingEdge, paste, collapse = ",")
    write.csv(fgsea_res, file = paste0("../GSEA_results/GSEA_", current, ".csv"), row.names = FALSE)
    
    # Plot top pathways
    topPathways <- head(fgsea_res, top_n)
    p <- ggplot(topPathways, aes(reorder(pathway, NES), NES)) +
      geom_col(aes(fill = padj < 0.05)) +
      coord_flip() + theme_minimal() +
      labs(title = current, y = "NES", x = "Pathway")
    ggsave(paste0("../GSEA_results/TopPathways_", current, ".png"), plot = p, width = 16, height = 12, dpi = 300)
  }
  
  # Print total number of successful results
  message("Final length of gsea_results: ", length(gsea_results))
  return(gsea_results)
}


gsea_results <- run_batch_fgsea(deg_results, pathways)

saveRDS(gsea_results, file = "../data/gsea_results.rds")

create_pathway_table_plot <- function(pathways, cluster, cont, n_top=10) {
  subset_name <- paste0(cluster, "_", gsub(" - ", "_", cont))
  fgseaRes <- gsea_results[[subset_name]]
  
  top10_DEGs <- read.csv("../DEG_results/Top10_DEGs.csv")
  ranks <- top10_DEGs$logFC
  names(ranks) <- top10_DEGs$gene
  
  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=n_top), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=n_top), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  table_plot <- plotGseaTable(pathways[topPathways], ranks, fgseaRes, 
                        gseaParam=0.5)
  ggsave(paste0("../GSEA_results/pathwayTable_", subset_name, ".png"), width = 16, height = 12, dpi = 300)
}

create_enrichment_plot <- function(pathways, target_pathway, n_top=10) {
  top10_DEGs <- read.csv("../DEG_results/Top10_DEGs.csv")
  ranks <- top10_DEGs$logFC
  names(ranks) <- top10_DEGs$gene
  
  plotEnrichment(pathways[[target_pathway]], ranks) + 
    labs(title="Programmed Cell Death",
         subtitle = target_pathway)
  
  ggsave(paste0("../GSEA_results/Enrichment_", target_pathway, ".png"))
}


# Plot pathway table
## args:

### required:
#### pathways: Pathways from database (has to be the same database for GSEA)
#### cluster: Cluster you are interested in (e.g. "B cell")
#### cont: Contrast (e.g. "PD - RD")

### optional:
#### n_top: Number of top genes you want to analyze

##Example: ("B cell" and "PD - RD")
create_pathway_table_plot(pathways = pathways, cluster = "B cell", cont = "PD - RD")


# Plot pathway table
## args:

### required:
#### pathways: Pathways from database (has to be the same database for GSEA)
#### target_pathway: Pathway you want to analyze

### optional:
#### n_top: Number of top genes you want to analyze

##Example: ("GOBP_ADAPTIVE_IMMUNE_RESPONSE")
create_enrichment_plot(pathways, target_pathway = "GOBP_ADAPTIVE_IMMUNE_RESPONSE")

