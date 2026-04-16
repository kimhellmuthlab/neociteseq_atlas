# ------------------------------------------------------------------------------
# Script: StackedBarPlot_cellproportions.R
# Input: Preprocessed single-cell dataset (.h5ad)
# Author: Paula Rothämel
# ------------------------------------------------------------------------------

# Description:
# This script loads CITE-seq data from an AnnData (.h5ad) object and generates
# stacked bar plots showing relative cell type proportions across gestational
# age groups. Reproduces Fig. 2B when applied to the myeloid cell subset.
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------------------------
library(reticulate)
library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(ggpattern)

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
# Prepare plot
# ------------------------------------------------------------------------------
meta <- seurat@meta.data

meta <- meta %>%
  dplyr::mutate(
    new_labels = case_when(
      cell_labels == "CD15+ myeloid cells" & 
        leiden_totalVI == "6" ~ "ARG1high CD15+",
      TRUE ~ cell_labels
    )
  )

prop_data <- meta %>%
  group_by(Age_Group, new_labels) %>%
  dplyr::summarise(count = n()) %>%
  group_by(Age_Group) %>%
  dplyr::mutate(prop = count / sum(count)) %>%
  ungroup()

# ------------------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------------------
fill_colors <- c(
  'Classical Monocytes' = '#4682B4', 
  'Intermediate Monocytes' = '#FF8C00', 
  'Non-classical Monocytes' = '#BA55D3',
  'cDCs' = '#A9A9A9',  
  'CD15+ myeloid cells' = '#9ACD32', 
  'ARG1high CD15+' = '#9ACD32'
)

patterns <- c(
  'Classical Monocytes' = 'none', 
  'Intermediate Monocytes' = 'none', 
  'Non-classical Monocytes' = 'none', 
  'CD15+ myeloid cells' = 'none', 
  "ARG1high CD15+" = 'stripe',
  'cDCs' = 'none')  

prop_data$new_labels <- factor(
  prop_data$new_labels,
  levels = c("ARG1high CD15+", "CD15+ myeloid cells", "Classical Monocytes", 
             "Intermediate Monocytes", "Non-classical Monocytes", "cDCs") 
)

ggplot(prop_data, aes(x = Age_Group, y = prop, fill = new_labels, 
                      pattern = new_labels)) +
  geom_bar_pattern(
    stat = "identity",
    color = "black",
    linewidth = 0.2,
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_density = 0.1,
    pattern_spacing = 0.04,
    width = 0.5
  ) +
  scale_fill_manual(values = fill_colors, name = "Cell Type") +
  scale_pattern_manual(values = patterns, guide = "none") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = c(0, 0)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        pattern = c("stripe", "none", "none", "none", "none", "none")
      )
    )
  ) +
  labs(y = "Proportion (%)", x = "") +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey80", linewidth = 0.3),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10)
  )

sessionInfo()
