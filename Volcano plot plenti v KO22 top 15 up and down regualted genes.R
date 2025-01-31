# Ensure gene names are set as rownames
rownames(results_plenti_vs_KO22) <- count_data$label
volcano_data_22 <- as.data.frame(results_plenti_vs_KO22)
volcano_data_22$Gene <- rownames(volcano_data_22)

# Define significance thresholds
volcano_data_22$Significance <- "Not Significant"
volcano_data_22$Significance[volcano_data_22$padj < 0.05 & volcano_data_22$log2FoldChange > log2(2)] <- "Upregulated"
volcano_data_22$Significance[volcano_data_22$padj < 0.05 & volcano_data_22$log2FoldChange < -log2(2)] <- "Downregulated"

# Separate top 15 upregulated and downregulated genes
top_upregulated_22 <- volcano_data_22 %>%
  filter(Significance == "Upregulated") %>%
  arrange(padj) %>%
  head(15)

top_downregulated_22 <- volcano_data_22 %>%
  filter(Significance == "Downregulated") %>%
  arrange(padj) %>%
  head(15)

# Combine top 30 genes
top_genes_22 <- rbind(top_upregulated_22, top_downregulated_22)

# Create the volcano plot
ggplot(volcano_data_22, aes(x = log2FoldChange, y = -log10(padj), color = Significance, label = Gene)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Not Significant" = "grey", "Upregulated" = "red", "Downregulated" = "blue")) + 
  geom_vline(xintercept = c(-log2(2), log2(2)), linetype = "dashed", color = "black") + #FC threshold
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + # p-value threshold
  geom_text(data = top_upregulated_22, aes(label = Gene), vjust = 1, hjust = 1, size = 3, color = "red", check_overlap = TRUE) + #Red for upregulated
  geom_text(data = top_downregulated_22, aes(label = Gene), vjust = 1, hjust = 1, size = 3, color = "blue", check_overlap = TRUE) + #Blue for downregulated
  xlim(-10, 10) +  # Limit the Log2 Fold Change range as it will go up to 20
  labs(title = "Volcano Plot (Plenti vs KO22) - Top 15 Up & Downregulated Genes", 
       x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme_minimal()
