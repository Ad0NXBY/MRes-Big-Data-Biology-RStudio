#PCA PLOTS
#Transform the data for PCA
vsd <- vst(dds, blind = FALSE) #variance stabilizing transformation

#Generate a PCA plot
plotPCA(vsd, intgroup = "condition") +
  ggtitle("PCA Plot of RNA-Seq Data")

#DESeq2 results for comparisons
#plenti vs KO22
results_plenti_vs_KO22 <- results(dds, contrast = c("condition", "KO22", "plenti"))
summary(results_plenti_vs_KO22) #summary of results
#plenti vs KO23
results_plenti_vs_KO23 <- results(dds, contrast = c("condition", "KO23", "plenti"))
summary(results_plenti_vs_KO23) #summary of results
