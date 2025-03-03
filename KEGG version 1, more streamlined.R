#### KEGG Enrichment Analysis ####
# Function to Convert Gene Symbols to Entrez IDs
convert_to_entrez <- function(gene_list) {
  entrez_ids <- mapIds(org.Hs.eg.db, 
                       keys = gene_list, 
                       column = "ENTREZID", 
                       keytype = "SYMBOL", 
                       multiVals = "first")
  
  # Remove NA values
  return(entrez_ids[!is.na(entrez_ids)])
}

# Function to Perform KEGG Enrichment
perform_kegg_enrichment <- function(gene_list) {
  entrez_ids <- convert_to_entrez(gene_list)  # Reuse the convert_to_entrez() function
  enrichKEGG(gene = entrez_ids,
             organism = "hsa",  # Human KEGG pathways
             pvalueCutoff = 0.05,
             keyType = "kegg")
}

# Function to Plot KEGG Results with ggplot2
plot_kegg_ggplot <- function(kegg_res, title, filename, top_n = 20) {
  if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res@result)) > 0) {
    kegg_df <- as.data.frame(kegg_res@result)
    kegg_df <- head(kegg_df[order(kegg_df$p.adjust), ], top_n)
    
    p <- ggplot(kegg_df, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
      geom_point() +
      scale_color_gradient(low = "red", high = "blue") +
      labs(title = title, x = "Gene Ratio", y = "KEGG Pathway", color = "Adjusted p-value", size = "Gene Count") +
      theme_minimal() +
      theme(text = element_text(size = 12))
    
    # Save plot
    ggsave(filename, plot = p, width = 10, height = 6)
  } else {
    print(paste("No significant KEGG pathways found for", title))
  }
}

# Convert Gene Symbols to Entrez IDs
up_entrez_KO22 <- convert_to_entrez(upregulated_genes_KO22)
down_entrez_KO22 <- convert_to_entrez(downregulated_genes_KO22)

# Perform KEGG Enrichment for Upregulated and Downregulated Genes
kegg_up_KO22 <- enrichKEGG(gene = up_entrez_KO22, organism = "hsa", pvalueCutoff = 0.05)
kegg_down_KO22 <- enrichKEGG(gene = down_entrez_KO22, organism = "hsa", pvalueCutoff = 0.05)

# Plot and Save KEGG Results
plot_kegg_ggplot(kegg_up_KO22, "KEGG Pathway for KO22 - Upregulated Genes", "KEGG_Upregulated_KO22.png")
plot_kegg_ggplot(kegg_down_KO22, "KEGG Pathway for KO22 - Downregulated Genes", "KEGG_Downregulated_KO22.png")

# Save Results
write.csv(as.data.frame(kegg_up_KO22@result), "KEGG_Upregulated_KO22.csv")
write.csv(as.data.frame(kegg_down_KO22@result), "KEGG_Downregulated_KO22.csv")

# Convert Gene Symbols to Entrez IDs
up_entrez_KO23 <- convert_to_entrez(upregulated_genes_KO23)
down_entrez_KO23 <- convert_to_entrez(downregulated_genes_KO23)

# Perform KEGG Enrichment for Upregulated and Downregulated Genes
kegg_up_KO23 <- enrichKEGG(gene = up_entrez_KO23, organism = "hsa", pvalueCutoff = 0.05)
kegg_down_KO23 <- enrichKEGG(gene = down_entrez_KO23, organism = "hsa", pvalueCutoff = 0.05)

# Plot and Save KEGG Results
plot_kegg_ggplot(kegg_up_KO23, "KEGG Pathway for KO23 - Upregulated Genes", "KEGG_Upregulated_KO23.png")
plot_kegg_ggplot(kegg_down_KO23, "KEGG Pathway for KO23 - Downregulated Genes", "KEGG_Downregulated_KO23.png")

# Save Results
write.csv(as.data.frame(kegg_up_KO23@result), "KEGG_Upregulated_KO23.csv")
write.csv(as.data.frame(kegg_down_KO23@result), "KEGG_Downregulated_KO23.csv")
