# Load necessary libraries for GO analysis for KO22
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)

# Define the function to convert fractions to decimals
fraction_to_decimal <- function(fraction_string) {
  sapply(strsplit(fraction_string, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
}

# Perform GO enrichment analysis for upregulated genes in different ontologies in KO22
upregulated_genes_KO22 <- rownames(results_plenti_vs_KO22)[results_plenti_vs_KO22$padj < 0.05 & results_plenti_vs_KO22$log2FoldChange > log2(2)]

enrich.go.up.bp_KO22 <- enrichGO(gene = upregulated_genes_KO22, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
enrich.go.up.cc_KO22 <- enrichGO(gene = upregulated_genes_KO22, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "CC", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
enrich.go.up.mf_KO22 <- enrichGO(gene = upregulated_genes_KO22, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "MF", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)

result.enrich.go.up.bp_KO22 <- enrich.go.up.bp_KO22@result
result.enrich.go.up.cc_KO22 <- enrich.go.up.cc_KO22@result
result.enrich.go.up.mf_KO22 <- enrich.go.up.mf_KO22@result

# Perform GO enrichment analysis for downregulated genes in different ontologies in KO22
downregulated_genes_KO22 <- rownames(results_plenti_vs_KO22)[results_plenti_vs_KO22$padj < 0.05 & results_plenti_vs_KO22$log2FoldChange < -log2(2)]

enrich.go.down.bp_KO22 <- enrichGO(gene = downregulated_genes_KO22, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
enrich.go.down.cc_KO22 <- enrichGO(gene = downregulated_genes_KO22, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "CC", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
enrich.go.down.mf_KO22 <- enrichGO(gene = downregulated_genes_KO22, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "MF", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)

result.enrich.go.down.bp_KO22 <- enrich.go.down.bp_KO22@result
result.enrich.go.down.cc_KO22 <- enrich.go.down.cc_KO22@result
result.enrich.go.down.mf_KO22 <- enrich.go.down.mf_KO22@result

# Convert GeneRatio from fraction to decimal for KO22
result.enrich.go.up.bp_KO22$GeneRatio <- fraction_to_decimal(result.enrich.go.up.bp_KO22$GeneRatio)
result.enrich.go.up.cc_KO22$GeneRatio <- fraction_to_decimal(result.enrich.go.up.cc_KO22$GeneRatio)
result.enrich.go.up.mf_KO22$GeneRatio <- fraction_to_decimal(result.enrich.go.up.mf_KO22$GeneRatio)

result.enrich.go.down.bp_KO22$GeneRatio <- fraction_to_decimal(result.enrich.go.down.bp_KO22$GeneRatio)
result.enrich.go.down.cc_KO22$GeneRatio <- fraction_to_decimal(result.enrich.go.down.cc_KO22$GeneRatio)
result.enrich.go.down.mf_KO22$GeneRatio <- fraction_to_decimal(result.enrich.go.down.mf_KO22$GeneRatio)

# Save GO enrichment results for KO22
write.csv(result.enrich.go.up.bp_KO22, "GO_enrichment_upregulated_BP_KO22.csv", row.names = FALSE)
write.csv(result.enrich.go.up.cc_KO22, "GO_enrichment_upregulated_CC_KO22.csv", row.names = FALSE)
write.csv(result.enrich.go.up.mf_KO22, "GO_enrichment_upregulated_MF_KO22.csv", row.names = FALSE)

write.csv(result.enrich.go.down.bp_KO22, "GO_enrichment_downregulated_BP_KO22.csv", row.names = FALSE)
write.csv(result.enrich.go.down.cc_KO22, "GO_enrichment_downregulated_CC_KO22.csv", row.names = FALSE)
write.csv(result.enrich.go.down.mf_KO22, "GO_enrichment_downregulated_MF_KO22.csv", row.names = FALSE)

# Filter and visualize top GO terms for KO22
top_15_up_bp_KO22 <- head(result.enrich.go.up.bp_KO22, 15)
top_15_up_cc_KO22 <- head(result.enrich.go.up.cc_KO22, 15)
top_15_up_mf_KO22 <- head(result.enrich.go.up.mf_KO22, 15)

top_15_down_bp_KO22 <- head(result.enrich.go.down.bp_KO22, 15)
top_15_down_cc_KO22 <- head(result.enrich.go.down.cc_KO22, 15)
top_15_down_mf_KO22 <- head(result.enrich.go.down.mf_KO22, 15)

# Dotplot for upregulated GO terms in BP for KO22
up_GO_dotplot_bp_KO22 <- ggplot(top_15_up_bp_KO22, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
  geom_point() +
  labs(title = "Top 15 Upregulated GO Terms (BP) KO22", x = "Gene Ratio", y = "GO Term") +
  theme_minimal()
ggsave("GO_upregulated_dotplot_BP_KO22.png", plot = up_GO_dotplot_bp_KO22, width = 10, height = 6)

# Dotplot for upregulated GO terms in CC for KO22
up_GO_dotplot_cc_KO22 <- ggplot(top_15_up_cc_KO22, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
  geom_point() +
  labs(title = "Top 15 Upregulated GO Terms (CC) KO22", x = "Gene Ratio", y = "GO Term") +
  theme_minimal()
ggsave("GO_upregulated_dotplot_CC_KO22.png", plot = up_GO_dotplot_cc_KO22, width = 10, height = 6)

# Dotplot for upregulated GO terms in MF for KO22
up_GO_dotplot_mf_KO22 <- ggplot(top_15_up_mf_KO22, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
  geom_point() +
  labs(title = "Top 15 Upregulated GO Terms (MF) KO22", x = "Gene Ratio", y = "GO Term") +
  theme_minimal()
ggsave("GO_upregulated_dotplot_MF_KO22.png", plot = up_GO_dotplot_mf_KO22, width = 10, height = 6)

# Dotplot for downregulated GO terms in BP for KO22
down_GO_dotplot_bp_KO22 <- ggplot(top_15_down_bp_KO22, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
  geom_point() +
  labs(title = "Top 15 Downregulated GO Terms (BP) KO22", x = "Gene Ratio", y = "GO Term") +
  theme_minimal()
ggsave("GO_downregulated_dotplot_BP_KO22.png", plot = down_GO_dotplot_bp_KO22, width = 10, height = 6)

# Dotplot for downregulated GO terms in CC for KO22
down_GO_dotplot_cc_KO22 <- ggplot(top_15_down_cc_KO22, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
  geom_point() +
  labs(title = "Top 15 Downregulated GO Terms (CC) KO22", x = "Gene Ratio", y = "GO Term") +
  theme_minimal()
ggsave("GO_downregulated_dotplot_CC_KO22.png", plot = down_GO_dotplot_cc_KO22, width = 10, height = 6)

# Dotplot for downregulated GO terms in MF for KO22
down_GO_dotplot_mf_KO22 <- ggplot(top_15_down_mf_KO22, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
  geom_point() +
  labs(title = "Top 15 Downregulated GO Terms (MF) KO22", x = "Gene Ratio", y = "GO Term") +
  theme_minimal()
ggsave("GO_downregulated_dotplot_MF_KO22.png", plot = down_GO_dotplot_mf_KO22, width = 10, height = 6)




# Perform GO enrichment analysis for upregulated genes in different ontologies in KO23
upregulated_genes_KO23 <- rownames(results_plenti_vs_KO23)[results_plenti_vs_KO23$padj < 0.05 & results_plenti_vs_KO23$log2FoldChange > log2(2)]

enrich.go.up.bp_KO23 <- enrichGO(gene = upregulated_genes_KO23, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
enrich.go.up.cc_KO23 <- enrichGO(gene = upregulated_genes_KO23, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "CC", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
enrich.go.up.mf_KO23 <- enrichGO(gene = upregulated_genes_KO23, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "MF", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)

result.enrich.go.up.bp_KO23 <- enrich.go.up.bp_KO23@result
result.enrich.go.up.cc_KO23 <- enrich.go.up.cc_KO23@result
result.enrich.go.up.mf_KO23 <- enrich.go.up.mf_KO23@result

# Perform GO enrichment analysis for downregulated genes in different ontologies in KO23
downregulated_genes_KO23 <- rownames(results_plenti_vs_KO23)[results_plenti_vs_KO23$padj < 0.05 & results_plenti_vs_KO23$log2FoldChange < -log2(2)]

enrich.go.down.bp_KO23 <- enrichGO(gene = downregulated_genes_KO23, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
enrich.go.down.cc_KO23 <- enrichGO(gene = downregulated_genes_KO23, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "CC", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
enrich.go.down.mf_KO23 <- enrichGO(gene = downregulated_genes_KO23, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "MF", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)

result.enrich.go.down.bp_KO23 <- enrich.go.down.bp_KO23@result
result.enrich.go.down.cc_KO23 <- enrich.go.down.cc_KO23@result
result.enrich.go.down.mf_KO23 <- enrich.go.down.mf_KO23@result

# Convert GeneRatio from fraction to decimal for KO23
result.enrich.go.up.bp_KO23$GeneRatio <- fraction_to_decimal(result.enrich.go.up.bp_KO23$GeneRatio)
result.enrich.go.up.cc_KO23$GeneRatio <- fraction_to_decimal(result.enrich.go.up.cc_KO23$GeneRatio)
result.enrich.go.up.mf_KO23$GeneRatio <- fraction_to_decimal(result.enrich.go.up.mf_KO23$GeneRatio)

result.enrich.go.down.bp_KO23$GeneRatio <- fraction_to_decimal(result.enrich.go.down.bp_KO23$GeneRatio)
result.enrich.go.down.cc_KO23$GeneRatio <- fraction_to_decimal(result.enrich.go.down.cc_KO23$GeneRatio)
result.enrich.go.down.mf_KO23$GeneRatio <- fraction_to_decimal(result.enrich.go.down.mf_KO23$GeneRatio)

# Save GO enrichment results for KO23
write.csv(result.enrich.go.up.bp_KO23, "GO_enrichment_upregulated_BP_KO23.csv", row.names = FALSE)
write.csv(result.enrich.go.up.cc_KO23, "GO_enrichment_upregulated_CC_KO23.csv", row.names = FALSE)
write.csv(result.enrich.go.up.mf_KO23, "GO_enrichment_upregulated_MF_KO23.csv", row.names = FALSE)

write.csv(result.enrich.go.down.bp_KO23, "GO_enrichment_downregulated_BP_KO23.csv", row.names = FALSE)
write.csv(result.enrich.go.down.cc_KO23, "GO_enrichment_downregulated_CC_KO23.csv", row.names = FALSE)
write.csv(result.enrich.go.down.mf_KO23, "GO_enrichment_downregulated_MF_KO23.csv", row.names = FALSE)

# Filter and visualize top GO terms for KO23
top_15_up_bp_KO23 <- head(result.enrich.go.up.bp_KO23, 15)
top_15_up_cc_KO23 <- head(result.enrich.go.up.cc_KO23, 15)
top_15_up_mf_KO23 <- head(result.enrich.go.up.mf_KO23, 15)

top_15_down_bp_KO23 <- head(result.enrich.go.down.bp_KO23, 15)
top_15_down_cc_KO23 <- head(result.enrich.go.down.cc_KO23, 15)
top_15_down_mf_KO23 <- head(result.enrich.go.down.mf_KO23, 15)

# Dotplot for upregulated GO terms in BP for KO23
up_GO_dotplot_bp_KO23 <- ggplot(top_15_up_bp_KO23, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
  geom_point() +
  labs(title = "Top 15 Upregulated GO Terms (BP) KO23", x = "Gene Ratio", y = "GO Term") +
  theme_minimal()
ggsave("GO_upregulated_dotplot_BP_KO23.png", plot = up_GO_dotplot_bp_KO23, width = 10, height = 6)

# Dotplot for upregulated GO terms in CC for KO23
up_GO_dotplot_cc_KO23 <- ggplot(top_15_up_cc_KO23, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
  geom_point() +
  labs(title = "Top 15 Upregulated GO Terms (CC) KO23", x = "Gene Ratio", y = "GO Term") +
  theme_minimal()
ggsave("GO_upregulated_dotplot_CC_KO23.png", plot = up_GO_dotplot_cc_KO23, width = 10, height = 6)

# Dotplot for upregulated GO terms in MF for KO23
up_GO_dotplot_mf_KO23 <- ggplot(top_15_up_mf_KO23, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
  geom_point() +
  labs(title = "Top 15 Upregulated GO Terms (MF) KO23", x = "Gene Ratio", y = "GO Term") +
  theme_minimal()
ggsave("GO_upregulated_dotplot_MF_KO23.png", plot = up_GO_dotplot_mf_KO23, width = 10, height = 6)

# Dotplot for downregulated GO terms in BP for KO23
down_GO_dotplot_bp_KO23 <- ggplot(top_15_down_bp_KO23, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
  geom_point() +
  labs(title = "Top 15 Downregulated GO Terms (BP) KO23", x = "Gene Ratio", y = "GO Term") +
  theme_minimal()
ggsave("GO_downregulated_dotplot_BP_KO23.png", plot = down_GO_dotplot_bp_KO23, width = 10, height = 6)

# Dotplot for downregulated GO terms in CC for KO23
down_GO_dotplot_cc_KO23 <- ggplot(top_15_down_cc_KO23, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
  geom_point() +
  labs(title = "Top 15 Downregulated GO Terms (CC) KO23", x = "Gene Ratio", y = "GO Term") +
  theme_minimal()
ggsave("GO_downregulated_dotplot_CC_KO23.png", plot = down_GO_dotplot_cc_KO23, width = 10, height = 6)

# Dotplot for downregulated GO terms in MF for KO23
down_GO_dotplot_mf_KO23 <- ggplot(top_15_down_mf_KO23, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
  geom_point() +
  labs(title = "Top 15 Downregulated GO Terms (MF) KO23", x = "Gene Ratio", y = "GO Term") +
  theme_minimal()
ggsave("GO_downregulated_dotplot_MF_KO23.png", plot = down_GO_dotplot_mf_KO23, width = 10, height = 6)

