# Install ggVennDiagram package if not installed
if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
  install.packages("ggVennDiagram")
}

# Load the package
library(ggVennDiagram)

DEG_KO22 <- results_plenti_vs_KO22$GeneSymbol
DEG_KO23 <- results_plenti_vs_KO23$GeneSymbol

# Create the Venn diagram
venn_data <- list(KO22_vs_plenti = DEG_KO22, KO23_vs_plenti = DEG_KO23)

# Plot the Venn diagram
ggVennDiagram(venn_data) +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal()
