# Ensure gene names are set as rownames
rownames(results_plenti_vs_KO23) <- count_data$label
volcano_data_23 <- as.data.frame(results_plenti_vs_KO23)
volcano_data_23$Gene <- rownames(volcano_data_23)

# Define significance thresholds
volcano_data_23$Significance <- "Not Significant"
volcano_data_23$Significance[volcano_data_23$padj < 0.05 & volcano_data_23$log2FoldChange > log2(2)] <- "Upregulated"
volcano_data_23$Significance[volcano_data_23$padj < 0.05 & volcano_data_23$log2FoldChange < -log2(2)] <- "Downregulated"

# Top 15 up/downregulated
top_upregulated_23 <- volcano_data_23 %>%
  filter(Significance == "Upregulated") %>%
  arrange(padj) %>%
  head(15)

top_downregulated_23 <- volcano_data_23 %>%
  filter(Significance == "Downregulated") %>%
  arrange(padj) %>%
  head(15)

# Combine top 30 genes
top_genes_23 <- rbind(top_upregulated_23, top_downregulated_23)

# Volcano plot for KO23
ggplot(volcano_data_23, aes(x = log2FoldChange, y = -log10(padj), color = Significance, label = Gene)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Not Significant" = "grey", "Upregulated" = "red", "Downregulated" = "blue")) + 
  geom_vline(xintercept = c(-log2(2), log2(2)), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text(data = top_upregulated_23, aes(label = Gene), vjust = 1, hjust = 1, size = 3, color = "red", check_overlap = TRUE) +
  geom_text(data = top_downregulated_23, aes(label = Gene), vjust = 1, hjust = 1, size = 3, color = "blue", check_overlap = TRUE) +
  labs(title = "Volcano Plot (Plenti vs KO23) - Top 15 Up & Downregulated Genes", 
       x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme_minimal()
