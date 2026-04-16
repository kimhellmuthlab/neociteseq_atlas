# ------------------------------------------------------------------------------
# Script: Pheatmap.R
# Input: Preprocessed single-cell dataset (.h5ad)
# Output: heatmap_myeloid_markers.pdf
# Author: Paula Rothämel
# ------------------------------------------------------------------------------

# Description:
# This script loads CITE-seq data from an AnnData (.h5ad) object, subsets CD15+
# myeloid cells, aggregates normalized RNA expression at the age-group level,
# scales gene expression, and visualizes selected marker genes as a heatmap.
# Samples are grouped by gestational age categories.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------------------------
library(reticulate)
library(Seurat)
library(pheatmap)
library(Matrix)
library(dplyr)

# ------------------------------------------------------------------------------
# Input / Output
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
input_file <- ifelse(length(args) >= 1, args[1], "PATH/TO/YOUR/FILE.h5ad")
output_file <- ifelse(length(args) >= 2, args[2], "heatmap_myeloid_markers.pdf")

gtf_file <- ifelse(
  length(args) >= 3,
  args[3],
  "data/gencode.v42.primary_assembly.annotation-filtered.gtf"
)

if (!file.exists(input_file)) stop("Input file not found: ", input_file)
if (!file.exists(gtf_file)) stop("GTF file not found: ", gtf_file)

# ------------------------------------------------------------------------------
# Helper Functions
# ------------------------------------------------------------------------------
extract_gtf_field <- function(x, field) {
  pattern <- paste0(field, " ([^;]+);")
  m <- regexpr(pattern, x)
  regmatches(x, m) |> sub(paste0(field, " "), "", x = _)
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
# Group by gestational age category
# ------------------------------------------------------------------------------
seurat$Age_Group <- ifelse(
  seurat$Sample_Name %in% c("unstim_late", "unstim_term"), ">=32w",
  ifelse(seurat$Sample_Name %in% c("unstim_extreme", "unstim_very"), "<32w", 
         "adult")
) 

# ------------------------------------------------------------------------------
# Map Ensembl IDs → gene symbols
# ------------------------------------------------------------------------------
gtf <- read.delim(
  gtf_file,
  header = FALSE,
  sep = "\t",
  comment.char = "#",
  stringsAsFactors = FALSE
)

gtf_gene <- gtf[gtf$V3 == "gene", ]

gene_id <- extract_field(gtf_gene$V9, "gene_id")
gene_name <- extract_field(gtf_gene$V9, "gene_name")
gene_id <- sub("\\..*$", "", gene_id)
gene_name <- sub(";$", "", gene_name)
ensg_to_symbol <- setNames(gene_name, gene_id)
counts <- GetAssayData(seurat, assay = "RNA", slot = "counts")
ensembl_ids <- rownames(counts)
ensembl_ids_clean <- sub("\\..*$", "", ensembl_ids)

gene_symbols <- ensg_to_symbol[ensembl_ids_clean]

# fallback for missing symbols
missing <- is.na(gene_symbols) | gene_symbols == ""
gene_symbols[missing] <- ensembl_ids_clean[missing]

rownames(counts) <- gene_symbols
counts <- counts[!duplicated(rownames(counts)), ]

# rebuild Seurat object
seurat <- CreateSeuratObject(
  counts = counts,
  meta.data = seurat@meta.data
)

seurat <- NormalizeData(seurat)

# ------------------------------------------------------------------------------
# Subset CD15+ myeloid cells
# ------------------------------------------------------------------------------
subset_obj <- subset(
  seurat,
  subset = cell_labels == "CD15+ myeloid cells"
)

subset_obj <- NormalizeData(subset_obj)

# ------------------------------------------------------------------------------
# Marker genes
# ------------------------------------------------------------------------------
markers.to.plot <- c("TNFAIP6", "CLEC4E", 
                     "CXCL1", "G0S2", "CXCL2", "IL1B", "ARG1", 
                     "ALPL", "FCGR3B", "CCL4L2", "CCL3L3", "CCL4", 
                     "IL1R2", "CCL20", "IL1A", "CXCL8", "CCL3")

# ------------------------------------------------------------------------------
# Heatmap
# ------------------------------------------------------------------------------
rna_data <- GetAssayData(subset_obj, assay = "RNA", slot = "data")[markers.to.plot,]
metadata <- subset_obj@meta.data

avg_expression <- sapply(markers.to.plot, function(gene) {
  tapply(rna_data[gene, ], metadata$Age_Group, mean, na.rm = TRUE)
})

avg_expression <- t(avg_expression)

avg_expression_scaled <- t(scale(t(avg_expression)))

pdf(output_file, width = 5, height = 5.5)

pheatmap(
  avg_expression_scaled, 
  cluster_rows = TRUE,  
  cluster_cols = FALSE, 
  scale = "none",  
  color = colorRampPalette(c("navy", "white", "#de2d26"))(100), 
  show_rownames = TRUE,  
  show_colnames = TRUE,
  angle_col = "45", 
  fontsize_row = 12.5,
  fontsize_col = 15 
)

dev.off()

sessionInfo()
