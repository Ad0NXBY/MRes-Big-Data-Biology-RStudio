# Ensure both data frames have Geneid as a column for matching
results_plenti_vs_KO22_df$Geneid <- rownames(results_plenti_vs_KO22_df)

# Merge to add gene labels to results_plenti_vs_KO22_df
results_plenti_vs_KO22_labeled <- merge(results_plenti_vs_KO22_df, count_data[, c("Geneid", "label")], 
                                        by = "Geneid", all.x = TRUE)

# Define significance thresholds
pval_cutoff <- 0.05
logfc_cutoff <- 2

# Update the significance column in the labeled data
results_plenti_vs_KO22_labeled$significance <- ifelse(
  results_plenti_vs_KO22_labeled$pvalue < pval_cutoff & abs(results_plenti_vs_KO22_labeled$log2FoldChange) > logfc_cutoff,
  "Significant",
  "Not Significant"
)

# Identify the most significant genes
most_significant <- results_plenti_vs_KO22_labeled %>%
  filter(pvalue < 0.001 & abs(log2FoldChange) > 2) %>%
  arrange(pvalue) %>%
  head(5)  # Select top 5 most significant genes

# Volcano plot with gene labels
ggplot(results_plenti_vs_KO22_labeled, aes(x = log2FoldChange, y = -log10(pvalue), color = significance)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("grey", "red")) +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "blue") +
  geom_text_repel(data = most_significant, aes(label = label), 
                  size = 3, max.overlaps = 10) +
  labs(title = "Volcano Plot of plenti vs KO22",
       x = "Log2 Fold Change",
       y = "-Log10 P-value") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
## Check chatgpt fixing label error to change it aroud
