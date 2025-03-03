#### KEGG Enrichment Analysis without Functions ####
# Convert Gene Symbols to Entrez IDs
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(DOSE)

# Convert Gene Symbols to Entrez IDs
up_entrez_KO22 <- mapIds(org.Hs.eg.db, keys = upregulated_genes_KO22, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
up_entrez_KO22 <- up_entrez_KO22[!is.na(up_entrez_KO22)]

down_entrez_KO22 <- mapIds(org.Hs.eg.db, keys = downregulated_genes_KO22, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
down_entrez_KO22 <- down_entrez_KO22[!is.na(down_entrez_KO22)]

# Perform KEGG Enrichment for Upregulated Genes
kegg_up_KO22 <- enrichKEGG(gene = up_entrez_KO22, organism = "hsa", pvalueCutoff = 0.05)
kegg_df_up_KO22 <- as.data.frame(kegg_up_KO22@result)
kegg_df_up_KO22 <- head(kegg_df_up_KO22[order(kegg_df_up_KO22$p.adjust), ], 20)

# Plot KEGG Results for Upregulated Genes
ggplot(kegg_df_up_KO22, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  labs(title = "KEGG Pathway for KO22 - Upregulated Genes", x = "Gene Ratio", y = "KEGG Pathway", color = "Adjusted p-value", size = "Gene Count") +
  theme_minimal() +
  theme(text = element_text(size = 12))
ggsave("KEGG_Upregulated_KO22.png", width = 10, height = 6)

# Perform KEGG Enrichment for Downregulated Genes
kegg_down_KO22 <- enrichKEGG(gene = down_entrez_KO22, organism = "hsa", pvalueCutoff = 0.05)
kegg_df_down_KO22 <- as.data.frame(kegg_down_KO22@result)
kegg_df_down_KO22 <- head(kegg_df_down_KO22[order(kegg_df_down_KO22$p.adjust), ], 20)

# Plot KEGG Results for Downregulated Genes
ggplot(kegg_df_down_KO22, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  labs(title = "KEGG Pathway for KO22 - Downregulated Genes", x = "Gene Ratio", y = "KEGG Pathway", color = "Adjusted p-value", size = "Gene Count") +
  theme_minimal() +
  theme(text = element_text(size = 12))
ggsave("KEGG_Downregulated_KO22.png", width = 10, height = 6)

# Save Results
write.csv(kegg_df_up_KO22, "KEGG_Upregulated_KO22.csv")
write.csv(kegg_df_down_KO22, "KEGG_Downregulated_KO22.csv")


# Convert Gene Symbols to Entrez IDs
up_entrez_KO23 <- mapIds(org.Hs.eg.db, keys = upregulated_genes_KO23, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
up_entrez_KO23 <- up_entrez_KO23[!is.na(up_entrez_KO23)]

down_entrez_KO23 <- mapIds(org.Hs.eg.db, keys = downregulated_genes_KO23, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
down_entrez_KO23 <- down_entrez_KO23[!is.na(down_entrez_KO23)]

# Perform KEGG Enrichment for Upregulated Genes
kegg_up_KO23 <- enrichKEGG(gene = up_entrez_KO23, organism = "hsa", pvalueCutoff = 0.05)
kegg_df_up_KO23 <- as.data.frame(kegg_up_KO23@result)
kegg_df_up_KO23 <- head(kegg_df_up_KO23[order(kegg_df_up_KO23$p.adjust), ], 20)

# Plot KEGG Results for Upregulated Genes
ggplot(kegg_df_up_KO23, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  labs(title = "KEGG Pathway for KO23 - Upregulated Genes", x = "Gene Ratio", y = "KEGG Pathway", color = "Adjusted p-value", size = "Gene Count") +
  theme_minimal() +
  theme(text = element_text(size = 12))
ggsave("KEGG_Upregulated_KO23.png", width = 10, height = 6)

# Perform KEGG Enrichment for Downregulated Genes
kegg_down_KO23 <- enrichKEGG(gene = down_entrez_KO23, organism = "hsa", pvalueCutoff = 0.05)
kegg_df_down_KO23 <- as.data.frame(kegg_down_KO23@result)
kegg_df_down_KO23 <- head(kegg_df_down_KO23[order(kegg_df_down_KO23$p.adjust), ], 20)

# Plot KEGG Results for Downregulated Genes
ggplot(kegg_df_down_KO23, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  labs(title = "KEGG Pathway for KO23 - Downregulated Genes", x = "Gene Ratio", y = "KEGG Pathway", color = "Adjusted p-value", size = "Gene Count") +
  theme_minimal() +
  theme(text = element_text(size = 12))
ggsave("KEGG_Downregulated_KO23.png", width = 10, height = 6)

# Save Results
write.csv(kegg_df_up_KO23, "KEGG_Upregulated_KO23.csv")
write.csv(kegg_df_down_KO23, "KEGG_Downregulated_KO23.csv")