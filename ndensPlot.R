# Load required libraries
library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(MuDataSeurat)
library(hdf5r)

# Load and normalize preprocessed data
seurat <- ReadH5MU("myeloid_baseline.h5mu")
DefaultAssay(seurat) <- "rna"
seurat <- NormalizeData(seurat)

# Assign age groups 
seurat$Age_Group <- ifelse(
  seurat$Sample_Name %in% c("unstim_late", "unstim_term"), ">=32w",
  ifelse(seurat$Sample_Name %in% c("unstim_extreme", "unstim_very"), "<32w", 
         "adult")
) 

# UMAP plot 
cell_type_colors <- c(
  'Classical Monocytes' = '#4682B4', # Yellow-Green  
  'Intermediate Monocytes' = '#FF8C00', # Dark Orange  
  'Non-classical Monocytes' = '#BA55D3', # Medium Orchid  
  'cDCs' = '#A9A9A9',  
  'CD15+ myeloid cells' = '#9ACD32'  # Steel Blue  
)

DimPlot(seurat, reduction = "UMAP", label = FALSE, group.by = "cell_labels") +
  scale_color_manual(values = cell_type_colors) + 
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 20),
    legend.key.size = unit(0.3, "cm"), 
  )

# ndens UMAP plots
group_conditions <- list(
  extreme_very = c("unstim_very", "unstim_extreme"), 
  late_term = c("unstim_late", "unstim_term"), 
  adult = "adult"
  )

# Create the list to store the density plots
plots_grouped_ndensity <- list()

Idents(seurat) <- "cell_labels"

# Loop through each group and create the density plots
for (group_name in names(group_conditions)) {
  # Get the conditions for the current group
  conditions <- group_conditions[[group_name]]
  # Subset the Seurat object based on the Sample_Name being one of the specified conditions
  subset <- subset(seurat, subset = Sample_Name %in% conditions)
  # Extract UMAP coordinates for the subset
  umap_coords <- as.data.frame(Embeddings(subset, reduction = "UMAP"))
  umap_coords$cell_type <- Idents(subset)  # Use renamed identities here
  umap_coords$Group <- group_name
  # Create the UMAP plot with normalized density contours (ndensity)
  p_ndensity <- ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = cell_type), alpha = 0.5) +  # Map color to cell_type
    stat_density_2d(aes(fill = after_stat(ndensity)), geom = "raster", contour = FALSE, alpha = 0.7) +  # Normalized density contours
    scale_color_manual(values = cell_type_colors) + 
    scale_fill_viridis_c() +
    theme(
      axis.text = element_blank(),  # Remove axis labels (numbers)
      axis.title = element_blank(),  # Remove axis titles
      axis.ticks = element_blank(),  # Remove axis ticks
      axis.line = element_line(color = "black"),  # Keep axis lines visible
      legend.position = "right",  # Move legend to the right
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Title styling
      legend.text = element_text(size = 10),  # Make legend text smaller
      legend.title = element_text(size = 20),  # Make legend title smaller
      legend.key.size = unit(0.5, "cm")  # Make legend keys smaller
    ) +
    ggtitle(paste("Normalized Density-UMAP", group_name))  # Add title with group name
  # Store the normalized density plot
  plots_grouped_ndensity[[group_name]] <- p_ndensity
}

# Display a specific normalized density plot (e.g., "adult")
plots_grouped_ndensity[["adult"]]
