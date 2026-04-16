# ------------------------------------------------------------------------------
# Script: DEGs_Analysis_Visualization.R
# Input: Preprocessed single-cell dataset (.h5ad)
# Author: Paula Rothämel
# ------------------------------------------------------------------------------

# Description:
# Analysis of CITE-seq single-cell RNA data (.h5ad) using Seurat and DESeq2.
# Includes pseudobulk differential expression, age-stratified and stimulus-specific 
# comparisons, enrichment analysis, interaction modeling, and visualization.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------------------------
library(reticulate)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(enrichR)
library(forcats)
library(tibble)
library(DESeq2)

# ------------------------------------------------------------------------------
# Input / Output
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
input_file <- ifelse(length(args) >= 1, args[1], "PATH/TO/YOUR/FILE.h5ad")

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

gene_id <- extract_gtf_field(gtf_gene$V9, "gene_id")
gene_name <- extract_gtf_field(gtf_gene$V9, "gene_name")
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


# Assign broader cell label categories within metadata 
seurat@meta.data$author_cell_labels <- 
  as.character(seurat@meta.data$author_cell_labels)
seurat@meta.data$author_cell_labels[seurat@meta.data$author_cell_labels %in% 
                                      c("Classical Monocytes", 
                                        "Non-classical Monocytes", 
                                        "Intermediate Monocytes")] <- "Monocytes"
seurat@meta.data$author_cell_labels <- factor(seurat@meta.data$author_cell_labels)

# Age-stratified DEG analysis e.g., for Myeloid cells (stim vs. baseline) Fig. 4
# -------------------------------
# SETTINGS
# -------------------------------
seurat$Sample_Donor <- paste(seurat$Sample_Name, seurat$donor_id, sep = "_")
seurat$new_batch <- ifelse(grepl("^adult_", seurat$Sample_Donor), "adult",
                           ifelse(grepl("_cbmc_", seurat$Sample_Donor) & 
                                    grepl("^unstim_", seurat$Sample_Donor), "cbmc_unstim",
                                  ifelse(grepl("_cbmc_", seurat$Sample_Donor) &
                                           grepl("^LPS_", seurat$Sample_Donor), "cbmc_LPS",
                                         ifelse(grepl("_cbmc_", seurat$Sample_Donor) & 
                                                  grepl("^R848_", seurat$Sample_Donor), "cbmc_R848",
                                                ifelse(grepl("_neo_|_NeoI|_NeoII", seurat$Sample_Donor) & 
                                                         grepl("^unstim_", seurat$Sample_Donor), "neo_unstim",
                                                       ifelse(grepl("_neo_|_NeoI|_NeoII", seurat$Sample_Donor) & 
                                                                grepl("^LPS_", seurat$Sample_Donor), "neo_LPS",
                                                              ifelse(grepl("_neo_|_NeoI|_NeoII", seurat$Sample_Donor) & 
                                                                       grepl("^R848_", seurat$Sample_Donor), "neo_R848",
                                                                     NA)))))))
  
cell_types <- c(
  "Monocytes",
  "CD15+ myeloid cells"
)

stimulus <- "LPS"
# stimulus <- "R848" # for R848 vs. baseline 
outdir <- "DESeq2_age_stratified_results"
age_groups <- c("lt32+0", "ge32+0")
dir.create(outdir, showWarnings = FALSE)

# -------------------------------
# SUBSET 
# -------------------------------

subset_all <- subset(
  seurat,
  subset = Sample_Name %in% c(
    "unstim_extreme", "unstim_very",
    "unstim_late", "unstim_term",
    paste0(stimulus, "_extreme"),
    paste0(stimulus, "_very"),
    paste0(stimulus, "_late"),
    paste0(stimulus, "_term")
  )
)

Idents(subset_all) <- "author_cell_labels"

# -------------------------------
# MAIN LOOP
# -------------------------------

for (cluster in cell_types) {
  
  message("Processing cell type: ", cluster)
  
  cluster_seurat <- subset(subset_all, author_cell_labels == cluster)
  
  if (ncol(cluster_seurat) < 30) {
    message("  -> too few cells, skipping")
    next
  }
  
  # -------------------------------
  # PSEUDOBULK
  # -------------------------------
  
  gene_counts <- AggregateExpression(
    cluster_seurat,
    group.by = c("new_batch", "Sample_Donor"),
    assays = "RNA",
    slot = "counts",
    return.seurat = FALSE
  )$RNA%>%
    as.data.frame() %>%
    rownames_to_column("gene")
  
  if (nrow(gene_counts) == 0) next
  
  counts_mat <- gene_counts %>% column_to_rownames("gene")
  
  # -------------------------------
  # BUILD colData 
  # -------------------------------
  
  sample_names <- colnames(counts_mat)
  
  colData_all <- data.frame(
    sample    = sample_names,
    condition = ifelse(grepl(stimulus, sample_names), stimulus, "baseline"),
    age_group = ifelse(grepl("extreme|very", sample_names), "lt32+0", "ge32+0"),
    batch     = ifelse(grepl("cbmc", sample_names), "cbmc", "pbmc")
  )
  
  rownames(colData_all) <- colData_all$sample
  
  colData_all$condition <- factor(colData_all$condition, levels = c("baseline", stimulus))
  colData_all$age_group <- factor(colData_all$age_group, levels = c("lt32+0", "ge32+0"))
  colData_all$batch     <- factor(colData_all$batch, levels = c("cbmc", "pbmc"))
  
  # -------------------------------
  # LOOP OVER AGE GROUPS
  # -------------------------------
  
  for (ag in age_groups) {
    
    message("  Age group: ", ag)
    
    keep_samples <- rownames(colData_all)[colData_all$age_group == ag]
    
    colData <- droplevels(colData_all[keep_samples, , drop = FALSE])
    countData <- counts_mat[, keep_samples, drop = FALSE]
    
    # need both baseline + stimulus
    if (length(unique(colData$condition)) < 2) {
      message("    -> only one condition present, skipping")
      next
    }
    
    # -------------------------------
    # DESEQ2
    # -------------------------------
    
    dds <- DESeqDataSetFromMatrix(
      countData = countData,
      colData   = colData,
      design    = ~ batch + condition
    )
    
    dds <- DESeq(dds)
    
    res <- results(
      dds,
      contrast = c("condition", stimulus, "baseline")
    )
    
    # -------------------------------
    # SAVE
    # -------------------------------
    
    prefix <- paste0(
      outdir, "/",
      gsub(" ", "_", cluster), "_",
      stimulus, "_",
      ag
    )
    
    saveRDS(dds, paste0(prefix, "_dds.rds"))
    saveRDS(res, paste0(prefix, "_stim_vs_baseline.rds"))
  }
}

# Enrichment dotplots (fig. S4)
dbs <- c("MSigDB_Hallmark_2020")

get_up_genes <- function(res, label = "") {
  df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    filter(!is.na(padj)) %>%
    filter(padj < 0.05, log2FoldChange > 1)  %>%
    filter(!grepl("^(ENSG|LINC)", gene))
  
  message(label, " → upregulated genes: ", nrow(df))
  df$gene
}

run_enrichr <- function(gene_list, cell, stimulus, age_group) {
  
  if (length(gene_list) < 5) return(NULL)
  
  enr <- enrichr(gene_list, dbs)[[1]]
  
  enr <- enr %>%
    filter(Adjusted.P.value < 0.05) %>%
    arrange(Adjusted.P.value)
  
  if (nrow(enr) > 10) {
    enr <- enr %>% slice_head(n = 10)
  }
  
  enr %>%
    mutate(
      log_padj = -log10(Adjusted.P.value),
      cell = cell,
      stimulus = stimulus,
      age_group = age_group
    )
}

build_enrichment_df <- function(
    res_lt32, res_ge32,
    cell, stimulus
) {
  
  message("Processing ", stimulus, " – ", cell)
  
  genes_lt32 <- get_up_genes(res_lt32, paste(cell, stimulus, "<32+0"))
  genes_ge32 <- get_up_genes(res_ge32, paste(cell, stimulus, "≥32+0"))
  
  bind_rows(
    run_enrichr(genes_lt32, cell, stimulus, "<32+0"),
    run_enrichr(genes_ge32, cell, stimulus, "≥32+0")
  )
}

plot_df <- bind_rows(
  
  ## R848 – Monocytes
  build_enrichment_df(
    readRDS("DESeq2_age_stratified_results/Monocytes_R848_lt32+0_stim_vs_baseline.rds"),
    readRDS("DESeq2_age_stratified_results/Monocytes_R848_ge32+0_stim_vs_baseline.rds"),
    "Monocytes", "R848"
  ),
  
  ## R848 – CD15+
  build_enrichment_df(
    readRDS("DESeq2_age_stratified_results/CD15+_myeloid_cells_R848_lt32+0_stim_vs_baseline.rds"),
    readRDS("DESeq2_age_stratified_results/CD15+_myeloid_cells_R848_ge32+0_stim_vs_baseline.rds"),
    "CD15+ myeloid cells", "R848"
  ),
  
  ## LPS – Monocytes
  build_enrichment_df(
    readRDS("DESeq2_age_stratified_results/Monocytes_LPS_lt32+0_stim_vs_baseline.rds"),
    readRDS("DESeq2_age_stratified_results/Monocytes_LPS_ge32+0_stim_vs_baseline.rds"),
    "Monocytes", "LPS"
  ),
  
  ## LPS – CD15+
  build_enrichment_df(
    readRDS("DESeq2_age_stratified_results/CD15+_myeloid_cells_LPS_lt32+0_stim_vs_baseline.rds"),
    readRDS("DESeq2_age_stratified_results/CD15+_myeloid_cells_LPS_ge32+0_stim_vs_baseline.rds"),
    "CD15+ myeloid cells", "LPS"
  )
)

plot_df$cell <- factor(
  plot_df$cell,
  levels = c("Monocytes", "CD15+ myeloid cells")
)

ggplot(
  plot_df,
  aes(
    x = age_group,
    y = fct_reorder(Term, log_padj),
    size = log_padj,
    color = age_group
  )
) +
  geom_point(alpha = 0.85) +
  facet_grid(stimulus ~ cell, scales = "free_y") +
  scale_size_continuous(range = c(2, 7), 
                        breaks = scales::pretty_breaks(n = 3)) +
  scale_color_manual(
    values = c("<32+0" = "#e16462", "≥32+0" = "#9a4fb5"),
    name = "Gestational age"
  ) +
  labs(
    x = "",
    y = "MSigDB Hallmark term",
    size = "-log10(p-adj)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )

# Prepare interaction analysis, scatter, and bar plots from age-stratified DESeq2-results

# -------------------------------
# LOAD AGE-STRATIFIED DESEQ2 RESULTS
# -------------------------------

load_res <- function(cell, stim, age){
  
  readRDS(
    paste0(
      "DESeq2_age_stratified_results/",
      cell, "_", stim, "_", age,
      "_stim_vs_baseline.rds"
    )
  )
}

# Monocytes
res_Mono_R848_lt32 <- load_res("Monocytes","R848","lt32+0")
res_Mono_R848_ge32 <- load_res("Monocytes","R848","ge32+0")

res_Mono_LPS_lt32 <- load_res("Monocytes","LPS","lt32+0")
res_Mono_LPS_ge32 <- load_res("Monocytes","LPS","ge32+0")

# CD15
res_CD15_R848_lt32 <- load_res("CD15+_myeloid_cells","R848","lt32+0")
res_CD15_R848_ge32 <- load_res("CD15+_myeloid_cells","R848","ge32+0")

res_CD15_LPS_lt32 <- load_res("CD15+_myeloid_cells","LPS","lt32+0")
res_CD15_LPS_ge32 <- load_res("CD15+_myeloid_cells","LPS","ge32+0")

# -------------------------------
# BUILD UNION DEG GENE SETS
# -------------------------------

make_union_geneset_from_res <- function(res_early, res_late){
  
  df_early <- as.data.frame(res_early) %>%
    tibble::rownames_to_column("gene")
  
  df_late <- as.data.frame(res_late) %>%
    tibble::rownames_to_column("gene")
  
  early_genes <- df_early %>%
    filter(!is.na(padj)) %>%
    filter(padj < 0.05, log2FoldChange > 1) %>%
    filter(!grepl("^ENSG|^LINC", gene)) %>%
    pull(gene)
  
  late_genes <- df_late %>%
    filter(!is.na(padj)) %>%
    filter(padj < 0.05, log2FoldChange > 1) %>%
    filter(!grepl("^ENSG|^LINC", gene)) %>%
    pull(gene)
  
  unique(c(early_genes, late_genes))
}

# -------------------------------
# CREATE GENE SETS FOR INTERACTION ANALYSIS
# -------------------------------

genes_monocytes_R848 <- make_union_geneset_from_res(
  res_Mono_R848_lt32,
  res_Mono_R848_ge32
)

genes_monocytes_LPS <- make_union_geneset_from_res(
  res_Mono_LPS_lt32,
  res_Mono_LPS_ge32
)

genes_cd15_R848 <- make_union_geneset_from_res(
  res_CD15_R848_lt32,
  res_CD15_R848_ge32
)

genes_cd15_LPS <- make_union_geneset_from_res(
  res_CD15_LPS_lt32,
  res_CD15_LPS_ge32
)

# Interaction DESeq2 analysis 
cell_types <- c(
  "Monocytes",
  "CD15+ myeloid cells"
)

stimulus <- "LPS"
# stimulus <- "R848"

outdir <- "DESeq2_interaction_results_genesets"
dir.create(outdir, showWarnings = FALSE)

# -------------------------------
# SUBSET ONCE: KEEP BOTH AGE GROUPS
# -------------------------------

subset_all <- subset(
  seurat,
  subset = Sample_Name %in% c(
    "unstim_extreme", "unstim_very",
    "unstim_late", "unstim_term",
    paste0(stimulus, "_extreme"),
    paste0(stimulus, "_very"),
    paste0(stimulus, "_late"),
    paste0(stimulus, "_term")
  )
)

Idents(subset_all) <- "author_cell_labels"

get_geneset <- function(cluster, stimulus) {
  if (cluster == "Monocytes" && stimulus == "R848") return(genes_monocytes_R848)
  if (cluster == "Monocytes" && stimulus == "LPS")  return(genes_monocytes_LPS)
  if (cluster == "CD15+ myeloid cells" && stimulus == "R848") return(genes_cd15_R848)
  if (cluster == "CD15+ myeloid cells" && stimulus == "LPS")  return(genes_cd15_LPS)
  NULL
}

# -------------------------------
# LOOP OVER CELL TYPES
# -------------------------------

for (cluster in cell_types) {
  
  message("Processing: ", cluster)
  
  cluster_seurat <- subset(subset_all, author_cell_labels == cluster)
  
  if (ncol(cluster_seurat) < 30) {
    message("  -> too few cells, skipping")
    next
  }
  
  # -------------------------------
  # PSEUDOBULK
  # -------------------------------
  
  gene_counts <- AggregateExpression(
    cluster_seurat,
    group.by = c("new_batch", "Sample_Donor"),
    assays = "RNA",
    slot = "counts",
    return.seurat = FALSE
  )$RNA %>%
    as.data.frame() %>%
    rownames_to_column("gene")
  
  if (nrow(gene_counts) == 0) next
  
  # -------------------------------
  # BUILD colData
  # -------------------------------
  
  sample_names <- colnames(gene_counts)[-1]
  
  condition <- ifelse(
    grepl(stimulus, sample_names),
    stimulus,
    "baseline"
  )
  
  age_group <- ifelse(
    grepl("extreme|very", sample_names),
    "lt32+0",
    "ge32+0"
  )
  
  batch <- ifelse(
    grepl("cbmc", sample_names),
    "cbmc",
    "pbmc"
  )
  
  colData <- data.frame(
    sample = sample_names,
    condition = condition,
    age_group = age_group,
    batch = batch
  )
  
  rownames(colData) <- colData$sample
  
  colData$condition <- factor(colData$condition, levels = c("baseline", stimulus))
  colData$age_group <- factor(colData$age_group, levels = c("lt32+0", "ge32+0"))
  colData$batch <- factor(colData$batch, levels = c("cbmc", "pbmc"))
  
  # ensure both conditions present
  if (length(unique(colData$condition)) < 2) {
    message("  -> only one condition present, skipping")
    next
  }
  
  # -------------------------------
  # RESTRICT TO GENE SET  <<< NEW
  # -------------------------------
  
  counts_mat <- gene_counts %>%
    column_to_rownames("gene") %>%
    as.matrix()
  
  geneset <- get_geneset(cluster, stimulus)
  
  if (is.null(geneset)) {
    message("  -> no gene set defined, skipping")
    next
  }
  
  counts_mat_gs <- counts_mat[
    rownames(counts_mat) %in% geneset,
    ,
    drop = FALSE
  ]
  
  if (nrow(counts_mat_gs) < 5) {
    message("  -> too few genes in gene set, skipping")
    next
  }
  
  # -------------------------------
  # DESEQ2 INTERACTION MODEL
  # -------------------------------
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts_mat_gs,
    colData = colData,
    design = ~ batch + age_group + condition + age_group:condition
  )
  
  dds <- DESeq(dds)
  
  # -------------------------------
  # RESULTS
  # -------------------------------
  
  # 1) Interaction term = STATISTICS
  res_interaction <- results(
    dds,
    name = paste0("age_groupge32.0.condition", stimulus)
  )
  
  # 2) Log2FC <32+0 (stim vs baseline)
  res_lt32 <- results(
    dds,
    contrast = c("condition", stimulus, "baseline")
  )
  
  # 3) Log2FC ≥32+0 (stim vs baseline)
  res_ge32 <- results(
    dds,
    list(
      c(
        paste0("condition_", stimulus, "_vs_baseline"),
        paste0("age_groupge32.0.condition", stimulus)
      )
    )
  )
  
  # -------------------------------
  # SAVE EVERYTHING
  # -------------------------------
  
  prefix <- paste0(outdir, "/", gsub(" ", "_", cluster), "_", stimulus)
  
  saveRDS(dds, paste0(prefix, "_dds.rds"))
  saveRDS(res_interaction, paste0(prefix, "_interaction.rds"))
  saveRDS(res_lt32, paste0(prefix, "_lt32_log2FC.rds"))
  saveRDS(res_ge32, paste0(prefix, "_ge32_log2FC.rds"))
}

# DEG scatter plot function 
prep_deseq_df <- function(res) {
  as.data.frame(res) %>%
    tibble::rownames_to_column("gene")
}

is_protein_coding <- function(gene) {
  !grepl("^ENSG", gene) & !grepl("^LINC", gene)
}

scatter_deg_union_realFC_interaction <- function(
    df_early,
    df_late,
    df_interaction,
    label_early,
    label_late,
    padj_cut = 0.05,
    log2fc_cut = 1,
    interaction_padj_cut = 0.1,
    top_n_interaction = 10
    
) {
  
  ## ---- significance (union rule)
  sig_early <- df_early %>%
    transmute(gene, sig_early = padj < padj_cut & log2FoldChange > log2fc_cut)
  
  sig_late <- df_late %>%
    transmute(gene, sig_late = padj < padj_cut & log2FoldChange > log2fc_cut)
  
  ## ---- FCs
  fc_early <- df_early %>%
    dplyr::select(gene, log2FC_early = log2FoldChange)
  
  fc_late <- df_late %>%
    dplyr::select(gene, log2FC_late = log2FoldChange)
  
  ## ---- interaction
  interaction <- df_interaction %>%
    dplyr::select(
      gene,
      interaction_padj = padj,
      interaction_log2FC = log2FoldChange
    )
  
  ## ---- merge
  merged <- full_join(fc_early, fc_late, by = "gene") %>%
    left_join(sig_early, by = "gene") %>%
    left_join(sig_late,  by = "gene") %>%
    left_join(interaction, by = "gene") %>%
    mutate(
      sig_early = replace_na(sig_early, FALSE),
      sig_late  = replace_na(sig_late,  FALSE),
      sig_category = sig_early | sig_late
    )
  
  ## ---- protein-coding only + inclusion
  merged <- merged %>%
    filter(is_protein_coding(gene)) %>%
    filter(sig_category)
  
  ## ---- interaction-based coloring
  merged <- merged %>%
    mutate(
      interaction_group = case_when(
        interaction_padj < interaction_padj_cut & interaction_log2FC > 0 ~ "higher in ≥32+0",
        interaction_padj < interaction_padj_cut & interaction_log2FC < 0 ~ "higher in <32+0",
        TRUE ~ "no interaction"
      )
    )
  
  ## ---- handle FC NAs
  merged <- merged %>%
    mutate(
      log2FC_early = replace_na(log2FC_early, 0),
      log2FC_late  = replace_na(log2FC_late,  0)
    )
  
  ## ---- LABELS: top N per direction, padj-ranked
  label_ge32 <- merged %>%
    filter(interaction_padj < interaction_padj_cut,
           interaction_log2FC > 0) %>%
    arrange(interaction_padj) %>%
    slice_head(n = top_n_interaction)
  
  label_lt32 <- merged %>%
    filter(interaction_padj < interaction_padj_cut,
           interaction_log2FC < 0) %>%
    arrange(interaction_padj) %>%
    slice_head(n = top_n_interaction)
  
  label_df <- bind_rows(label_ge32, label_lt32)
  
  ## ---- plot
  p <- ggplot(
    merged,
    aes(
      x = log2FC_late,
      y = log2FC_early,
      color = interaction_group, 
      alpha = interaction_group
    )
  ) +
    geom_point(size = 0.7) +
    geom_text_repel(
      data = label_df,
      aes(label = gene),
      size = 3,
      show.legend = FALSE, 
      max.overlaps = 8
    ) +
    geom_abline(slope = 1, linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    scale_color_manual(
      values = c(
        "higher in <32+0" = "#e16462",
        "higher in ≥32+0" = "#9a4fb5",
        "no interaction" = "grey75"
      )
    ) +
    scale_alpha_manual(
      values = c(
        "higher in <32+0" = 1,
        "higher in ≥32+0" = 1,
        "no interaction" = 0.5
      ),
      guide = "none"   # Alpha nicht in der Legende anzeigen
    ) +
    labs(
      x = paste0("log2FC of ", label_late),
      y = paste0("log2FC of ", label_early),
      color = "Stronger induction"
    ) +
    theme_minimal(base_size = 12)
  
  list(
    plot  = p,
    genes = merged$gene,
    table = merged
  )
}

# e.g., for Monocytes R848 vs. baseline 
res_Mono_R848_lt32 <- readRDS(
  "DESeq2_age_stratified_results/Monocytes_R848_lt32+0_stim_vs_baseline.rds")
res_Mono_R848_ge32 <- readRDS(
  "DESeq2_age_stratified_results/Monocytes_R848_ge32+0_stim_vs_baseline.rds")
res_Mono_R848_interaction <- readRDS(
  "DESeq2_interaction_results_genesets/Monocytes_R848_interaction.rds")

out_Mono_R848 <- scatter_deg_union_realFC_interaction(
  prep_deseq_df(res_Mono_R848_lt32),
  prep_deseq_df(res_Mono_R848_ge32),
  prep_deseq_df(res_Mono_R848_interaction),
  "<32+0", "≥32+0"
)

p_Mono_R848     <- out_Mono_R848$plot
genes_Mono_R848 <- out_Mono_R848$genes
tbl_Mono_R848   <- out_Mono_R848$table

p_Mono_R848

# for Monocytes LPS vs. baseline 
res_Mono_LPS_lt32 <- readRDS(
  "DESeq2_age_stratified_results/Monocytes_LPS_lt32+0_stim_vs_baseline.rds")

res_Mono_LPS_ge32 <- readRDS(
  "DESeq2_age_stratified_results/Monocytes_LPS_ge32+0_stim_vs_baseline.rds")

res_Mono_LPS_interaction <- readRDS(
  "DESeq2_interaction_results_genesets/Monocytes_LPS_interaction.rds")


out_Mono_LPS <- scatter_deg_union_realFC_interaction(
  prep_deseq_df(res_Mono_LPS_lt32),
  prep_deseq_df(res_Mono_LPS_ge32),
  prep_deseq_df(res_Mono_LPS_interaction),
  "<32+0", "≥32+0"
)

p_Mono_LPS     <- out_Mono_LPS$plot
genes_Mono_LPS <- out_Mono_LPS$genes
tbl_Mono_LPS   <- out_Mono_LPS$table

p_Mono_LPS

# for CD15 R848 vs. baseline 
res_CD15_R848_lt32 <- readRDS(
  "DESeq2_age_stratified_results/CD15+_myeloid_cells_R848_lt32+0_stim_vs_baseline.rds")

res_CD15_R848_ge32 <- readRDS(
  "DESeq2_age_stratified_results/CD15+_myeloid_cells_R848_ge32+0_stim_vs_baseline.rds")

res_CD15_R848_interaction <- readRDS(
  "DESeq2_interaction_results_genesets/CD15+_myeloid_cells_R848_interaction.rds")


out_CD15_R848 <- scatter_deg_union_realFC_interaction(
  prep_deseq_df(res_CD15_R848_lt32),
  prep_deseq_df(res_CD15_R848_ge32),
  prep_deseq_df(res_CD15_R848_interaction),
  "<32+0", "≥32+0"
)

p_CD15_R848     <- out_CD15_R848$plot
genes_CD15_R848 <- out_CD15_R848$genes
tbl_CD15_R848   <- out_CD15_R848$table

p_CD15_R848

# for CD15 LPS vs. baseline 
res_CD15_LPS_lt32 <- readRDS(
  "DESeq2_age_stratified_results/CD15+_myeloid_cells_LPS_lt32+0_stim_vs_baseline.rds")

res_CD15_LPS_ge32 <- readRDS(
  "DESeq2_age_stratified_results/CD15+_myeloid_cells_LPS_ge32+0_stim_vs_baseline.rds")

res_CD15_LPS_interaction <- readRDS(
  "DESeq2_interaction_results_genesets/CD15+_myeloid_cells_LPS_interaction.rds")


out_CD15_LPS <- scatter_deg_union_realFC_interaction(
  prep_deseq_df(res_CD15_LPS_lt32),
  prep_deseq_df(res_CD15_LPS_ge32),
  prep_deseq_df(res_CD15_LPS_interaction),
  "<32+0", "≥32+0"
)

p_CD15_LPS     <- out_CD15_LPS$plot
genes_CD15_LPS <- out_CD15_LPS$genes
tbl_CD15_LPS   <- out_CD15_LPS$table

p_CD15_LPS

# -------------------------------
# BUILD TABLE OF TOP INTERACTION GENES
# -------------------------------

final_scatter_table <- bind_rows(
  tbl_Mono_R848 %>% mutate(cell = "Monocytes", stimulus = "R848"),
  tbl_CD15_R848 %>% mutate(cell = "CD15+ myeloid cells", stimulus = "R848"),
  tbl_Mono_LPS  %>% mutate(cell = "Monocytes", stimulus = "LPS"),
  tbl_CD15_LPS  %>% mutate(cell = "CD15+ myeloid cells", stimulus = "LPS")
)

top5_interaction_table <- final_scatter_table %>%
  filter(!is.na(interaction_padj)) %>%
  filter(interaction_padj < 0.1) %>%
  mutate(
    direction = case_when(
      interaction_log2FC > 0 ~ "higher in ≥32+0",
      interaction_log2FC < 0 ~ "higher in <32+0"
    )
  ) %>%
  group_by(cell, stimulus, direction) %>%
  arrange(interaction_padj) %>%
  slice_head(n = 5) %>%
  ungroup()

# DEG barplots, separate for R848 and LPS
deseq_lookup <- list(
  
  Monocytes = list(
    R848 = list(
      lt32 = res_Mono_R848_lt32,
      ge32 = res_Mono_R848_ge32,
      interaction = res_Mono_R848_interaction
    ),
    LPS = list(
      lt32 = res_Mono_LPS_lt32,
      ge32 = res_Mono_LPS_ge32,
      interaction = res_Mono_LPS_interaction
    )
  ),
  
  `CD15+ myeloid cells` = list(
    R848 = list(
      lt32 = res_CD15_R848_lt32,
      ge32 = res_CD15_R848_ge32,
      interaction = res_CD15_R848_interaction
    ),
    LPS = list(
      lt32 = res_CD15_LPS_lt32,
      ge32 = res_CD15_LPS_ge32,
      interaction = res_CD15_LPS_interaction
    )
  )
)

get_genes_from_top5 <- function(cell, stimulus, direction) {
  top5_interaction_table %>%
    filter(
      cell == !!cell,
      stimulus == !!stimulus,
      direction == !!direction
    ) %>%
    arrange(interaction_padj) %>%
    pull(gene)
}

extract_fc_and_stats <- function(cell, stimulus, genes) {
  res_lt32 <- deseq_lookup[[cell]][[stimulus]]$lt32
  res_ge32 <- deseq_lookup[[cell]][[stimulus]]$ge32
  tibble(
    gene = genes,
    log2FC_lt32 = res_lt32[genes, "log2FoldChange"],
    lfcSE_lt32  = res_lt32[genes, "lfcSE"],
    log2FC_ge32 = res_ge32[genes, "log2FoldChange"],
    lfcSE_ge32  = res_ge32[genes, "lfcSE"]
  )
}

make_gene_barplot_top5_single_stimulus <- function(
    cell_type,
    stimulus_keep,
    direction_keep,
    outfile
) {
  genes_of_interest <- get_genes_from_top5(
    cell = cell_type,
    stimulus = stimulus_keep,
    direction = direction_keep
  )
  plot_df <- extract_fc_and_stats(
    cell     = cell_type,
    stimulus = stimulus_keep,
    genes    = genes_of_interest
  )
  
  plot_long <- plot_df %>%
    pivot_longer(
      cols = c(log2FC_lt32, log2FC_ge32, lfcSE_lt32, lfcSE_ge32),
      names_to = c(".value", "age_group"),
      names_pattern = "(log2FC|lfcSE)_(lt32|ge32)"
    ) %>%
    mutate(
      age_group = recode(age_group, lt32 = "<32+0", ge32 = "≥32+0"),
      stimulus  = factor(stimulus_keep, levels = stimulus_keep),
      gene      = factor(gene, levels = genes_of_interest),
      
      ci_lower = log2FC - qnorm(0.975) * lfcSE,
      ci_upper = log2FC + qnorm(0.975) * lfcSE
    )
  
  p <- ggplot(
    plot_long,
    aes(x = stimulus, y = log2FC, fill = age_group)
  ) +
    geom_col(
      position = position_dodge(width = 0.5),
      width = 0.6
    ) +
    geom_errorbar(
      aes(ymin = ci_lower, ymax = ci_upper),
      position = position_dodge(width = 0.5),
      width = 0.2,
      linewidth = 0.2
    ) +
    facet_wrap(~ gene, nrow = 1, scales = "fixed") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    scale_fill_manual(values = c("<32+0" = "#e16462", "≥32+0" = "#9a4fb5")) +
    labs(
      x = "",
      y = "log2FC (stimulus vs. baseline)",
      fill = "Gestational age"
    ) +
    theme_classic(base_size = 9) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 9),
      axis.line = element_line(color = "grey40", linewidth = 0.3),
      axis.ticks = element_line(color = "grey40", linewidth = 0.3),
      axis.text.x = element_text(size = 9),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 9)
    )
  
  ggsave(
    filename = outfile,
    plot = p,
    width = 3.3,
    height = 2,
    units = "in"
  )
  
  return(p)
}

# e.g., for Monocytes R848 vs baseline - higher in <32+0 
make_gene_barplot_top5_single_stimulus(
  cell_type = "Monocytes",
  stimulus_keep = "R848",
  direction_keep = "higher in <32+0",
  outfile = "Monocytes_R848_top5_lt32.pdf"
)

sessionInfo()
