# Load required libraries
library(Seurat)
library(pheatmap)
library(MuDataSeurat)
library(hdf5r)

# Load data, e.g., baseline Myeloid cells 
seurat <- ReadH5MU("myeloid_baseline.h5mu")

# Assign age groups 
seurat$Age_Group <- ifelse(
  seurat$Sample_Name %in% c("unstim_late", "unstim_term"), ">=32w",
  ifelse(seurat$Sample_Name %in% c("unstim_extreme", "unstim_very"), "<32w", 
         "adult")
) 

# Subset and normalize data 
subset <- subset(seurat, manual_labels %in% c("CD15+ myeloid cells"))
DefaultAssay(subset) <- "rna"
subset <- NormalizeData(subset)

# e.g. for heatmap Fig. 2C
markers.to.plot <- c("TNFAIP6", "CLEC4E", 
                     "CXCL1", "G0S2", "CXCL2", "IL1B", "ARG1", 
                     "ALPL", "FCGR3B", "CCL4L2", "CCL3L3", "CCL4", 
                     "IL1R2", "CCL20", "IL1A", "CXCL8", "CCL3")

# Extract normalized expression data for the selected genes
rna_data <- GetAssayData(subset, assay = "rna", layer = "data")[markers.to.plot,]

# Extract the metadata
metadata <- subset@meta.data

# Calculate the average expression for each gene per group
avg_expression <- sapply(markers.to.plot, function(gene) {
  # Calculate the average expression of each gene across the different groups
  tapply(rna_data[gene, ], metadata$Age_Group, mean, na.rm = TRUE)
})

# Transpose for better readability
avg_expression <- t(avg_expression)

# View the average expression
head(avg_expression)

# Scale the matrix by rows (genes)
avg_expression_scaled <- t(scale(t(avg_expression)))

# Pheatmap 
pheatmap(
  avg_expression_scaled, 
  cluster_rows = TRUE,  
  cluster_cols = FALSE, 
  scale = "none", # already scaled by rows (genes) 
  color = colorRampPalette(c("navy", "white", "#de2d26"))(100), 
  show_rownames = TRUE,  
  show_colnames = TRUE,
  angle_col = "45",        # Spaltenbeschriftung rotieren
  fontsize_row = 12.5,       # Schriftgröße Zeilen
  fontsize_col = 15 # Schriftgröße Spalten
)
