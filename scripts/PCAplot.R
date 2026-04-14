# ------------------------------------------------------------------------------
# Script: PCAplot.R
# Input: Preprocessed single-cell dataset (.h5ad)
# Output: PCA_plot.pdf
# Author: Paula Rothämel
# ------------------------------------------------------------------------------

# Description:
# This script loads CITE-seq data from an AnnData (.h5ad) object, aggregates RNA
# expression at the sample level, performs PCA, and visualizes samples annotated
# by gestational age, treatment condition, and postnatal age.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------------------------
library(reticulate)
library(Seurat)
library(Matrix)
library(matrixStats)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(svglite)
library(grid)
library(ggnewscale)

# ------------------------------------------------------------------------------
# Input / Output
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
input_file <- ifelse(length(args) >= 1, args[1], "PATH/TO/YOUR/FILE.h5ad")
output_file <- ifelse(length(args) >= 2, args[2], "PCA_plot.pdf")

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
# Add metadata
# ------------------------------------------------------------------------------
meta <- py_to_r(adata$obs)
meta <- as.data.frame(meta)
seurat <- AddMetaData(seurat, metadata = meta)

# ------------------------------------------------------------------------------
# Aggregate expression
# ------------------------------------------------------------------------------
seurat$Sample_Donor <- paste(seurat$Sample_Name, seurat$donor_id, sep = "_")
seurat <- NormalizeData(seurat)

sample_expression_RNA <- AggregateExpression(
  seurat,
  group.by = "Sample_Donor",
  assays = "RNA",
  layer = "data",
  return.seurat = FALSE
)

aggregated_matrix <- sample_expression_RNA$RNA

# Remove zero-variance genes
filtered_matrix <- aggregated_matrix[apply(aggregated_matrix, 1, var) > 0, ]

if (nrow(filtered_matrix) == 0) {
  stop("No genes with non-zero variance found")
}

# ------------------------------------------------------------------------------
# PCA
# ------------------------------------------------------------------------------
set.seed(123)
pca_result <- prcomp(t(as.matrix(filtered_matrix)), scale. = TRUE)
explained_variance <- summary(pca_result)$importance[2, ] * 100

pca_df <- as.data.frame(pca_result$x)
pca_df$Sample_Donor <- rownames(pca_df)

# ------------------------------------------------------------------------------
# Metadata parsing
# ------------------------------------------------------------------------------
pca_df$Sample_Donor <- gsub("-", "_", pca_df$Sample_Donor)

pca_df <- pca_df %>%
  mutate(
    donor_group = case_when(
      grepl("unstim_|LPS_|R848_", Sample_Donor) ~ sub(".*_(neo_donor\\d+|cbmc_donor\\d+)", "\\1", Sample_Donor),
      TRUE ~ Sample_Donor
    ),
    treatment_group = case_when(
      grepl("^LPS_", Sample_Donor) ~ "LPS",
      grepl("^R848_", Sample_Donor) ~ "R848",
      TRUE ~ "baseline"
    ),
    Gestational_Age = case_when(
      grepl("extreme", tolower(Sample_Donor)) ~ "extreme",
      grepl("very", tolower(Sample_Donor)) ~ "very",
      grepl("late", tolower(Sample_Donor)) ~ "late",
      grepl("term", tolower(Sample_Donor)) ~ "term",
      grepl("adult", tolower(Sample_Donor)) ~ "adult",
      TRUE ~ NA_character_
    )
  )

pca_df$cell_type <- ifelse(grepl("cbmc", pca_df$donor_group), "CBMC", "PBMC")

pca_df <- pca_df %>%
  mutate(
    Postnatal_Age = case_when(
      cell_type == "CBMC" ~ "birth",
      cell_type == "PBMC" & Gestational_Age %in% c("extreme", "very", "late", "term") ~ "DOL3",
      cell_type == "PBMC" & Gestational_Age == "adult" ~ "adult",
      TRUE ~ NA_character_
    )
  )


# ------------------------------------------------------------------------------
# Plot 
# ------------------------------------------------------------------------------

# Colors and shapes
custom_colors <- c(
  "adult" = "#0d0887",
  "term" = "#6a00a8",
  "late" = "#b12a90",
  "very" = "#e16462",
  "extreme" = "#fca636"
)

custom_shapes <- c(
  "birth" = 24,
  "DOL3" = 21,
  "adult" = 22
)

# Ensure factors
pca_df$Gestational_Age <- factor(
  pca_df$Gestational_Age,
  levels = c("extreme", "very", "late", "term", "adult")
)

pca_df$Postnatal_Age <- factor(
  pca_df$Postnatal_Age,
  levels = c("birth", "DOL3", "adult")
)

# Condition variable
pca_df$Condition <- ifelse(pca_df$treatment_group == "baseline", "Baseline", "Stimulated")

# Labels (R / L)
pca_df$label <- dplyr::case_when(
  pca_df$treatment_group == "R848" ~ "R",
  pca_df$treatment_group == "LPS"  ~ "L",
  TRUE ~ NA_character_
)

# ------------------------------------------------------------------------------
# Fix PCA orientation for reproducibility
# ------------------------------------------------------------------------------
# Note:
# Principal components are only defined up to a sign (i.e., PC axes can be
# mirrored without changing the underlying data structure). Small numerical
# differences (e.g., due to data import via AnnData/reticulate or matrix
# representations) may result in flipped PC directions across runs.
#
# We fix the orientation here to ensure consistent and reproducible visualization.
pca_df$PC2 <- -pca_df$PC2

# Split data
baseline_df <- dplyr::filter(pca_df, Condition == "Baseline")
stim_df     <- dplyr::filter(pca_df, Condition == "Stimulated")

# ------------------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------------------

p <- ggplot() +
  
  # -------------------------
# Baseline (hollow points)
# -------------------------
geom_point(
  data = baseline_df,
  aes(x = PC1, y = PC2, color = Gestational_Age, shape = Postnatal_Age),
  fill = "white",
  size = 3.5,
  stroke = 1
) +
  
  # -------------------------
# Stimulated (filled points)
# -------------------------
geom_point(
  data = stim_df,
  aes(x = PC1, y = PC2, fill = Gestational_Age, color = Gestational_Age, shape = Postnatal_Age),
  size = 3.5,
  stroke = 1
) +
  
  # Labels (L / R)
  geom_text(
    data = pca_df,
    aes(x = PC1, y = PC2, label = label),
    size = 2,
    fontface = "bold",
    na.rm = TRUE
  ) +
  
  # -------------------------
# Main scales
# -------------------------
scale_color_manual(
  values = custom_colors,
  name = "Gestational Age"
) +
  
  scale_fill_manual(
    values = custom_colors,
    guide = "none"
  ) +
  
  scale_shape_manual(
    values = custom_shapes,
    name = "Postnatal Age"
  ) +
  
  # -------------------------
# New fill scale for condition legend
# -------------------------
ggnewscale::new_scale_fill() +
  
  # Dummy points for condition legend
  geom_point(
    data = data.frame(
      Condition = factor(c("Baseline", "LPS", "R848"),
                         levels = c("Baseline", "LPS", "R848")),
      x = NA, y = NA
    ),
    aes(x = x, y = y, fill = Condition),
    shape = 21,
    size = 3.5,
    color = "black",
    stroke = 1,
    inherit.aes = FALSE
  ) +
  
  scale_fill_manual(
    name = "Condition",
    values = c(
      "Baseline" = "white",
      "LPS" = "black",
      "R848" = "black"
    )
  ) +
  
  labs(
    title = "RNA PCA Plot of Samples",
    x = paste0("PC1 (", round(explained_variance[1], 1), "%)"),
    y = paste0("PC2 (", round(explained_variance[2], 1), "%)")
  ) +
  
  theme_light()

p

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------
ggsave(output_file, plot = p, device = "pdf", width = 5, height = 4)
write.csv(pca_df, "PCA_coordinates.csv", row.names = FALSE)

message("PCA plot saved to: ", output_file)
sessionInfo()
