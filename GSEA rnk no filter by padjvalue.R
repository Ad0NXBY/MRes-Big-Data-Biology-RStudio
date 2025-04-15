results_plenti_vs_KO22 <- results(dds, contrast = c("condition", "KO22", "plenti")) %>%
  as.data.frame() %>%
  na.omit()  # Remove genes with NA in padj

# Map Ensembl IDs to gene symbols
results_plenti_vs_KO22$GeneSymbol <- mapIds(
  org.Hs.eg.db,
  keys = rownames(results_plenti_vs_KO22),
  column = "SYMBOL",
  keytype = "SYMBOL",
  multiVals = "first"
)

# Remove genes without valid symbols
results_plenti_vs_KO22_clean <- results_plenti_vs_KO22 %>% 
  dplyr::filter(!is.na(GeneSymbol))

# Keep the most significant entry per gene
results_plenti_vs_KO22 <- results_plenti_vs_KO22_clean %>%
  group_by(GeneSymbol) %>%
  arrange(padj, desc(abs(log2FoldChange))) %>%
  slice(1) %>%
  ungroup()

# Convert to data frame, then set gene symbols as row names
results_plenti_vs_KO22 <- as.data.frame(results_plenti_vs_KO22)
rownames(results_plenti_vs_KO22) <- results_plenti_vs_KO22$GeneSymbol



# Sort the DEG list by log2FoldChange (highest to lowest)
results_plenti_vs_KO22 <- results_plenti_vs_KO22 %>%
  dplyr::arrange(desc(log2FoldChange))  # Sort by log2FoldChange

# Save DEG list as a .csv file
write.csv(
  results_plenti_vs_KO22, 
  file.path(data_dir, "DEG_plenti_vs_KO22.csv"), 
  row.names = TRUE
)

# Prepare ranked list for GSEA
ranked_list_KO22 <- results_plenti_vs_KO22 %>%
  dplyr::arrange(desc(log2FoldChange)) %>%  # Ensure it's sorted by highest to lowest log2FoldChange
  dplyr::select(GeneSymbol, log2FoldChange)

# Convert to a data frame
ranked_list_KO22 <- as.data.frame(ranked_list_KO22)

# Set gene symbols as row names and remove the redundant column
rownames(ranked_list_KO22) <- ranked_list_KO22$GeneSymbol
ranked_list_KO22 <- ranked_list_KO22 %>% dplyr::select(-GeneSymbol)

# Write .rnk file for GSEA
write.table(
  ranked_list_KO22,
  file.path(data_dir, "GSEA_ranked_list_KO22_nofilter_by_padj_.rnk"),
  sep = "\t",
  col.names = FALSE,
  row.names = TRUE,
  quote = FALSE
)





results_plenti_vs_KO23 <- results(dds, contrast = c("condition", "KO23", "plenti")) %>%
  as.data.frame() %>%
  na.omit()  # Remove genes with NA in padj

# Map Ensembl IDs to gene symbols
results_plenti_vs_KO23$GeneSymbol <- mapIds(
  org.Hs.eg.db,
  keys = rownames(results_plenti_vs_KO23),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Remove genes without valid symbols
results_plenti_vs_KO23_clean <- results_plenti_vs_KO23 %>% 
  dplyr::filter(!is.na(GeneSymbol))

# Keep the most significant entry per gene
results_plenti_vs_KO23 <- results_plenti_vs_KO23_clean %>%
  group_by(GeneSymbol) %>%
  arrange(padj, desc(abs(log2FoldChange))) %>%
  slice(1) %>%
  ungroup()

# Convert to data frame, then set gene symbols as row names
results_plenti_vs_KO23 <- as.data.frame(results_plenti_vs_KO23)
rownames(results_plenti_vs_KO23) <- results_plenti_vs_KO23$GeneSymbol



# Sort the DEG list by log2FoldChange (highest to lowest)
results_plenti_vs_KO23 <- results_plenti_vs_KO23 %>%
  dplyr::arrange(desc(log2FoldChange))  # Sort by log2FoldChange

# Save DEG list as a .csv file
write.csv(
  results_plenti_vs_KO23, 
  file.path(data_dir, "DEG_plenti_vs_KO23.csv"), 
  row.names = TRUE
)

# Prepare ranked list for GSEA
ranked_list_KO23 <- results_plenti_vs_KO23 %>%
  dplyr::arrange(desc(log2FoldChange)) %>%  # Ensure it's sorted by highest to lowest log2FoldChange
  dplyr::select(GeneSymbol, log2FoldChange)

# Convert to a data frame
ranked_list_KO23 <- as.data.frame(ranked_list_KO23)

# Set gene symbols as row names and remove the redundant column
rownames(ranked_list_KO23) <- ranked_list_KO23$GeneSymbol
ranked_list_KO23 <- ranked_list_KO23 %>% dplyr::select(-GeneSymbol)

# Write .rnk file for GSEA
write.table(
  ranked_list_KO23,
  file.path(data_dir, "GSEA_ranked_list_KO23_nofilter_by_padj_.rnk"),
  sep = "\t",
  col.names = FALSE,
  row.names = TRUE,
  quote = FALSE
)

