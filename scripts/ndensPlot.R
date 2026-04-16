# ------------------------------------------------------------------------------
# Script: ndensPlot.R
# Input: Preprocessed single-cell dataset (.h5ad)
# Author: Paula Rothämel
# ------------------------------------------------------------------------------

# Description:
# This script loads CITE-seq data from an AnnData (.h5ad) object, extracts 
# pre-computed UMAP embeddings, subsets by gestational age groups, and generates 
# normalized density (ndensity) UMAP plots overlaid with cell type annotations.
# Reproduces Fig. 2A when applied to the myeloid cell subset.
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------------------------
library(reticulate)
library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)

# ------------------------------------------------------------------------------
# Input
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
input_file <- ifelse(length(args) >= 1, args[1], "PATH/TO/YOUR/FILE.h5ad")

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# ------------------------------------------------------------------------------
# Load .h5ad via anndata
# ------------------------------------------------------------------------------
anndata <- import("anndata", convert = FALSE)
adata <- anndata$read_h5ad(input_file)

# ------------------------------------------------------------------------------
# Build Seurat object
# ------------------------------------------------------------------------------
# --- counts from raw.X ---
counts <- py_to_r(adata$raw$X$toarray())   # dense Python CSR → R matrix
counts <- as(counts, "dgCMatrix")          # sparse dgCMatrix for Seurat
counts <- t(counts)                        # gene × cell

# Fix rownames/colnames
rownames(counts) <- py_to_r(adata$raw$var_names$to_list())
colnames(counts) <- py_to_r(adata$obs_names$to_list())
seurat <- CreateSeuratObject(counts = counts)

# ------------------------------------------------------------------------------
# Load UMAP embeddings from AnnData
# ------------------------------------------------------------------------------
umap <- py_to_r(adata$obsm[["X_umap"]])
umap <- as.matrix(umap)
rownames(umap) <- colnames(seurat)

# Add to Seurat 
seurat[["umap"]] <- CreateDimReducObject(
  embeddings = umap,
  key = "UMAP_",
  assay = DefaultAssay(seurat)
)

# ------------------------------------------------------------------------------
# Add metadata and normalize
# ------------------------------------------------------------------------------
meta <- py_to_r(adata$obs)
meta <- as.data.frame(meta)
seurat <- AddMetaData(seurat, metadata = meta)

seurat <- NormalizeData(seurat)

# ------------------------------------------------------------------------------
# Group by gestational age category
# ------------------------------------------------------------------------------
seurat$Age_Group <- ifelse(
  seurat$Sample_Name %in% c("unstim_late", "unstim_term"), ">=32w",
  ifelse(seurat$Sample_Name %in% c("unstim_extreme", "unstim_very"), "<32w", 
         "adult")
) 

# ------------------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------------------

cell_type_colors <- c(
  'Classical Monocytes' = '#4682B4', 
  'Intermediate Monocytes' = '#FF8C00', 
  'Non-classical Monocytes' = '#BA55D3',   
  'cDCs' = '#A9A9A9',  
  'CD15+ myeloid cells' = '#9ACD32'   
)

# UMAP plot 
DimPlot(seurat, reduction = "umap", label = FALSE, group.by = "author_cell_labels") +
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

Idents(seurat) <- "author_cell_labels"

# Loop through each group and create the density plots
for (group_name in names(group_conditions)) {
  # Get the conditions for the current group
  conditions <- group_conditions[[group_name]]
  # Subset the Seurat object based on the Sample_Name being one of the specified conditions
  seurat_subset <- subset(
    seurat,
    subset = !is.na(Sample_Name) & Sample_Name %in% conditions
  )
  # Extract UMAP coordinates for the subset
  umap_coords <- as.data.frame(Embeddings(seurat_subset, reduction = "umap"))
  umap_coords$cell_type <- Idents(seurat_subset)
  umap_coords$Group <- group_name
  # Create the UMAP plot with normalized density contours (ndensity)
  p_ndensity <- ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = cell_type), alpha = 0.5) +  
    stat_density_2d(aes(fill = after_stat(ndensity)), geom = "raster", contour = FALSE, alpha = 0.7) +  
    scale_color_manual(values = cell_type_colors) + 
    scale_fill_viridis_c() +
    theme(
      axis.text = element_blank(), 
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_line(color = "black"),
      legend.position = "right",
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 20),
      legend.key.size = unit(0.5, "cm")
    ) +
    ggtitle(paste("Normalized Density-UMAP", group_name))
  # Store the normalized density plot
  plots_grouped_ndensity[[group_name]] <- p_ndensity
}

# Display a specific normalized density plot (e.g., "adult")
plots_grouped_ndensity[["adult"]]

sessionInfo()
