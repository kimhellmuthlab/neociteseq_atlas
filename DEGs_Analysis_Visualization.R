# Load required libraries
library(Seurat)
library(Matrix)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(MuDataSeurat)
library(hdf5r)

# Load data
# seurat <- ReadH5MU("mudata_baseline_all.h5mu")

# Assign broader cell label categories 


# Merge seurat objects


# Adapt metadata 
merged_obj$Stage <- ifelse(
  merged_obj$Sample_Name == "adult", 
  "adult", 
  "neonatal"
)

merged_obj$Stage_Celltype <- paste(merged_obj$Stage, merged_obj$cell_labels, sep = "_")

# Generate Heatmap 
# Normalize data 
merged_obj <- NormalizeData(merged_obj, verbose = FALSE)

# For e.g. Fig 2L 
markers.to.plot <- consistently_down_filtered

# Extract normalized expression data for the selected genes
rna_data <- GetAssayData(merged_obj, assay = "rna", layer = "data")[markers.to.plot,]

# Extract the metadata
metadata <- merged_obj@meta.data

# Calculate the average expression for each gene per group
avg_expression <- sapply(markers.to.plot, function(gene) {
  # Calculate the average expression of each gene across the different groups
  tapply(rna_data[gene, ], metadata$Stage_Celltype, mean, na.rm = TRUE)
}) # The result is a matrix with the average expression values for each gene and group

# Transpose for better readability
avg_expression <- t(avg_expression)

# View the average expression
head(avg_expression)

# Scale the matrix by rows (genes)
avg_expression_scaled <- t(scale(t(avg_expression)))

# Reorder rows based on the desired population order if you wanna show the groups in a special order
desired_population_order <- c(
  "neonatal_Monocytes",
  "adult_Monocytes", 
  "neonatal_CD15+",
  "adult_CD15+", 
  "neonatal_NK cells",
  "adult_NK cells", 
  "neonatal_CD4+ T cells",
  "adult_CD4+ T cells",
  "neonatal_CD8+ T cells",
  "adult_CD8+ T cells", 
  "neonatal_B cells",
  "adult_B cells" 
)

# Reorder the columns based on the population order
avg_expression_scaled <- avg_expression_scaled [, desired_population_order]

pheatmap(
  avg_expression_scaled, 
  cluster_rows = TRUE,  
  cluster_cols = FALSE, 
  scale = "none", # data already scaled 
  color = colorRampPalette(c("navy", "white", "#de2d26"))(100), 
  show_rownames = FALSE,  
  show_colnames = TRUE,
  angle_col = "45", 
  fontsize_row = 12.5,      
  fontsize_col = 15 
)