# ------------------------------------------------------------------------------
# Script: PCAplot.R
# Input: Preprocessed single-cell dataset (e.g., CELLxGENE dataset “All - neonatal and adult circulating immune cells baseline and stim”)  
# Output: PCA_plot_all.pdf
# Author: Paula Rothämel
# ------------------------------------------------------------------------------

# Description:
# This script loads CITE-seq data from an AnnData (.h5ad)object, aggregates RNA expression
# at the sample level, performs PCA, and visualizes samples annotated by
# gestational age, treatment condition, and postnatal age.
# ------------------------------------------------------------------------------

# Load required libraries
library(Seurat)
library(SeuratDisk)
library(Matrix)
library(matrixStats)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(svglite)
library(grid)

# Requires manual download of the CellxGene dataset “All - neonatal and adult circulating immune cells baseline and stim”
# The dataset should be placed in the working directory before running this script.
# This script aggregates RNA expression by sample, performs PCA, and generates the Fig. 1B plot.

input_file <- "mudata_all.h5mu"
output_file <- "PCA_plot_all.pdf"

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# Load and normalize data
seurat <- ReadH5MU(input_file)
DefaultAssay(seurat) <- "rna"
seurat <- NormalizeData(seurat)

seurat$Sample_Donor <- paste(seurat$Sample_Name, seurat$donor_id, sep = "_")

# Aggregate RNA data by donor and condition
sample_expression_RNA <- AggregateExpression(
  seurat,
  group.by = "Sample_Donor",    
  assays = "rna", 
  slot = "data", 
  return.seurat = FALSE     
)

# Extract the aggregated matrix
aggregated_matrix <- sample_expression_RNA$rna
filtered_matrix <- aggregated_matrix[apply(aggregated_matrix, 1, var) > 0, ]

# Perform PCA on the filtered data (samples as rows, genes as columns)
pca_result <- prcomp(t(filtered_matrix), scale. = TRUE)  # Transpose to have samples as rows and genes as columns

explained_variance <- summary(pca_result)$importance[2, ] * 100  # Get percentage

# Create a data frame for PCA coordinates
pca_df <- as.data.frame(pca_result$x)
pca_df$Sample_Donor <- rownames(pca_df)

# Standardize Sample_Name format by replacing hyphens with underscores in both pca_df and Seurat metadata
pca_df$Sample_Donor <- gsub("-", "_", pca_df$Sample_Donor)

# Create new grouping variables
pca_df <- pca_df %>%
  mutate(
    donor_group = case_when(
      grepl("unstim_|LPS_|R848_", Sample_Donor) ~ sub(".*_(neo_donor\\d+|cbmc_donor\\d+)", "\\1", Sample_Donor),
      TRUE ~ Sample_Donor
    ),
    treatment_group = case_when(
      grepl("^LPS_", Sample_Donor) ~ "LPS",
      grepl("^R848_", Sample_Donor) ~ "R848",
      TRUE ~ "unstim"
    )
  )

pca_df <- pca_df %>%
  mutate(
    Gestational_Age = case_when(
      grepl("extreme", Sample_Donor) ~ "extreme",
      grepl("very", Sample_Donor) ~ "very",
      grepl("late", Sample_Donor) ~ "late",
      grepl("term", Sample_Donor) ~ "term",
      grepl("adult", Sample_Donor) ~ "adult",
      TRUE ~ NA_character_  # Assign NA if no match is found
    )
  )

pca_df$cell_type <- ifelse(grepl("cbmc", pca_df$donor_group), "CBMC", "PBMC")

pca_df <- pca_df %>%
  mutate(
    Postnatal_Age = case_when(
      cell_type == "CBMC" ~ "birth",  # alle CBMC sind von Geburt
      cell_type == "PBMC" & Gestational_Age %in% c("extreme", "very", "late", "term") ~ "DOL3",
      cell_type == "PBMC" & Gestational_Age == "adult" ~ "adult",
      TRUE ~ NA_character_
    )
  )

# PCA plot 
custom_colors <- c(
  "adult"   = "#0d0887",
  "term"    = "#6a00a8",
  "late"    = "#b12a90",
  "very"    = "#e16462",
  "extreme" = "#fca636"
)

# Shapes for sample types
custom_shapes <- c(
  "birth" = 24,  # circle with border
  "DOL3" = 21, # triangle with border
  "adult" = 22
)

# Map fill: white for unstim, Age color for stim
pca_df$fill_color <- ifelse(
  pca_df$treatment_group == "unstim",
  "white",
  # z. B. grau für R848
  custom_colors[pca_df$Gestational_Age]  # für LPS
)

# Map outline color always to Age
pca_df$outline_color <- custom_colors[pca_df$Gestational_Age]

# Add a label column: "R" for R848, "L" for LPS, "" otherwise
pca_df$label <- ifelse(
  pca_df$treatment_group == "R848", "R",
  ifelse(pca_df$treatment_group == "LPS", "L", "")
)

p <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(
    aes(fill = fill_color, shape = Postnatal_Age, color = outline_color),
    size = 3.5, stroke = 0.8
  ) +
  geom_text(
    aes(label = label),
    vjust = 0.5, hjust = 0.5,
    size = 2, fontface = "bold"
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_shape_manual(values = custom_shapes) +
  guides(
    fill = guide_legend(
      title = "Age",
      override.aes = list(shape = 21, size = 3.5, stroke = 0.8)
    ),
    shape = guide_legend(
      title = "Postnatal Age",
      override.aes = list(fill = "gray70", color = "black", stroke = 0.8)
    )
  ) +
  labs(
    title = "RNA PCA Plot of Samples",
    x = paste0("PC1 (", round(explained_variance[1], 1), "%)"),
    y = paste0("PC2 (", round(explained_variance[2], 1), "%)")
  ) +
  theme_light() +
  theme(
    plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9),
    legend.key.size = unit(1.2, "lines"),
    legend.key.height = unit(1.5, "lines"),
    panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
    panel.grid.minor = element_blank()
  )

# save plot 
ggsave(filename = output_file, plot = p, device = "pdf", width = 5, height = 4)
