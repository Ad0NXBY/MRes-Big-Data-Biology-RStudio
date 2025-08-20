library(GSVA)
library(msigdbr)
library(GSEABase)

BiocManager::install("GSEABase")


# Get human gene sets from msigdbr
msigdbr_species <- msigdbr(species = "Homo sapiens")
hallmark_genesets_df <- msigdbr_species[msigdbr_species$gs_cat == "H", ]

# Convert to GeneSetCollection
library(GSEABase)

# Split gene sets
gset.idx.list <- split(hallmark_genesets_df$gene_symbol, hallmark_genesets_df$gs_name)

# Create GeneSet objects with unique genes per set
geneSets <- lapply(names(gset.idx.list), function(name) {
  GeneSet(geneIds = unique(gset.idx.list[[name]]), 
          setName = name, 
          geneIdType = SymbolIdentifier())
})

geneSetCollection <- GeneSetCollection(geneSets)

expr_matrix <- assay(vsd)  # if you've done variance stabilizing transformation


# Create parameter object and run GSVA
param <- gsvaParam(expr = expr_matrix, 
                   geneSets = geneSetCollection, 
                   kcdf = "Poisson")
gsva_results <- gsva(param)


head(gsva_results)

write.csv(gsva_results, file.path(data_dir, "GSVA_results.csv"), row.names = TRUE)



# Calculate t-test statistics for Control vs KO22
t_results_KO22 <- apply(gsva_results, 1, function(row) {
  t.test(row[c("FIH22KO1_Count", "FIH22KO2_Count", "FIH22KO3_Count")], 
         row[c("Control1_Count", "Control2_Count", "Control3_Count")])$statistic
})

# Remove "HALLMARK_" prefix from gene set names
names(t_results_KO22) <- sub("HALLMARK_", "", names(t_results_KO22))

# Create data frame for plotting
bar_data_KO22 <- data.frame(
  GeneSet = names(t_results_KO22),
  T_Value = unlist(t_results_KO22),
  FillColor = ifelse(abs(unlist(t_results_KO22)) < 1, "grey", ifelse(unlist(t_results_KO22) > 0, "red", "blue"))
)

# Sort by T value
bar_data_KO22 <- bar_data_KO22[order(bar_data_KO22$T_Value, decreasing = TRUE),]

# Plot bar graph
ggplot(bar_data_KO22, aes(x = reorder(GeneSet, T_Value), y = T_Value, fill = FillColor)) +
  geom_bar(stat = "identity") +
  scale_fill_identity() +
  coord_flip() +
  labs(x = "Gene Set", y = "t value of GSVA", title = "GSVA Analysis for Control vs KO22") +
  theme_minimal() +
  theme(text = element_text(color = "black"),
        legend.position = "none")




# Calculate t-test statistics for Control vs KO23
t_results_KO23 <- apply(gsva_results, 1, function(row) {
  t.test(row[c("FIH23KO1_Count", "FIH23KO2_Count", "FIH23KO3_Count")], 
         row[c("Control1_Count", "Control2_Count", "Control3_Count")])$statistic
})

# Remove "HALLMARK_" prefix from gene set names
names(t_results_KO23) <- sub("HALLMARK_", "", names(t_results_KO23))

# Create data frame for plotting
bar_data_KO23 <- data.frame(
  GeneSet = names(t_results_KO23),
  T_Value = unlist(t_results_KO23),
  FillColor = ifelse(abs(unlist(t_results_KO23)) < 0.05, "grey", ifelse(unlist(t_results_KO23) > 0, "red", "blue"))
)

# Sort by T value
bar_data_KO23 <- bar_data_KO23[order(bar_data_KO23$T_Value, decreasing = TRUE),]

# Plot bar graph
ggplot(bar_data_KO23, aes(x = reorder(GeneSet, T_Value), y = T_Value, fill = FillColor)) +
  geom_bar(stat = "identity") +
  scale_fill_identity() +
  coord_flip() +
  labs(x = "Gene Set", y = "t value of GSVA", title = "GSVA Analysis for Control vs KO23") +
  theme_minimal() +
  theme(text = element_text(color = "black"),
        legend.position = "none")



# Calculate t-test statistics for KO22 vs KO23
t_results_ko22_vs_ko23 <- apply(gsva_results, 1, function(row) {
  t.test(row[c("FIH22KO1_Count", "FIH22KO2_Count", "FIH22KO3_Count")],
         row[c("FIH23KO1_Count", "FIH23KO2_Count", "FIH23KO3_Count")])$statistic
})

# Remove "HALLMARK_" prefix from gene set names
names(t_results_ko22_vs_ko23) <- sub("HALLMARK_", "", names(t_results_ko22_vs_ko23))

# Create data frame for plotting
bar_data_ko22_vs_ko23 <- data.frame(
  GeneSet = names(t_results_ko22_vs_ko23),
  T_Value = unlist(t_results_ko22_vs_ko23),
  FillColor = ifelse(abs(unlist(t_results_ko22_vs_ko23)) < 1, "grey",
                     ifelse(unlist(t_results_ko22_vs_ko23) > 0, "red", "blue"))
)

# Sort by T value
bar_data_ko22_vs_ko23 <- bar_data_ko22_vs_ko23[order(bar_data_ko22_vs_ko23$T_Value, decreasing = TRUE),]

# Plot the GSVA t-values for KO22 vs KO23
ggplot(bar_data_ko22_vs_ko23, aes(x = reorder(GeneSet, T_Value), y = T_Value, fill = FillColor)) +
  geom_bar(stat = "identity") +
  scale_fill_identity() +
  coord_flip() +
  labs(x = "Gene Set", y = "t value of GSVA", title = "GSVA: KO22 vs KO23") +
  theme_minimal() +
  theme(text = element_text(color = "black"),
        legend.position = "none")





# Extract relevant GSVA scores for all 3 conditions
gsva_subset <- gsva_results[, c("Control1_Count", "Control2_Count", "Control3_Count",
                                "FIH22KO1_Count", "FIH22KO2_Count", "FIH22KO3_Count",
                                "FIH23KO1_Count", "FIH23KO2_Count", "FIH23KO3_Count")]

# Optionally, remove HALLMARK_ prefix
rownames(gsva_subset) <- sub("HALLMARK_", "", rownames(gsva_subset))

# Optionally scale rows (mean-center each gene set)
gsva_scaled <- t(scale(t(gsva_subset)))

# Create annotation for columns
annotation_col <- data.frame(
  Condition = c(rep("Control", 3), rep("KO22", 3), rep("KO23", 3))
)
rownames(annotation_col) <- colnames(gsva_scaled)


library(ComplexHeatmap)

# Create HeatmapAnnotation for color coding the conditions
ha <- HeatmapAnnotation(
  Condition = annotation_col$Condition,
  col = list(Condition = c("Control" = "#4DBBD5", "KO22" = "#E64B35", "KO23" = "#00A087"))
)

# Create the heatmap using ComplexHeatmap
GSVA_heatmap <- Heatmap(gsva_scaled,  # GSVA results matrix
        name = "GSVA Score",  # Colorbar label
        top_annotation = ha,  # Annotations for columns (conditions)
        cluster_rows = TRUE,  # Cluster rows (gene sets)
        cluster_columns = FALSE,  # Don't cluster columns to preserve order (optional)
        show_row_names = TRUE,  # Show row (gene set) names
        show_column_names = TRUE,  # Show column (sample) names
        col = colorRampPalette(c("navy", "white", "firebrick3"))(100),  # Custom color palette
        heatmap_legend_param = list(title = "GSVA Score"),  # Color bar title
        row_names_gp = gpar(fontsize = 7)# Adjust row name font size
)

jpeg(file.path(plot_dir,"GSVA_heatmap.jpg"), width = 1200, height = 1000, res = 150)
draw(GSVA_heatmap)
dev.off()

