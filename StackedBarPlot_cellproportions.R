# Load required libraries
library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(MuDataSeurat)
library(hdf5r)
library(ggpattern)

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

# Stacked barplot to visualize cell proportions
meta <- seurat@meta.data

meta <- meta %>%
  dplyr::mutate(
    new_labels = case_when(
      cell_labels == "CD15+ myeloid cells" & leiden_totalVI == "6" ~ "ARG1high CD15+",
      TRUE ~ cell_labels
    )
  )

prop_data <- meta %>%
  group_by(Age_Group, new_labels) %>%
  dplyr::summarise(count = n()) %>%
  group_by(Age_Group) %>%
  dplyr::mutate(prop = count / sum(count)) %>%  # sum(count) over same Age_Group only
  ungroup()

fill_colors <- c(
  'Classical Monocytes' = '#4682B4', # Yellow-Green  
  'Intermediate Monocytes' = '#FF8C00', # Dark Orange  
  'Non-classical Monocytes' = '#BA55D3', # Medium Orchid  
  'cDCs' = '#A9A9A9',  
  'CD15+ myeloid cells' = '#9ACD32', 
  'ARG1high CD15+' = '#9ACD32'# Steel Blue  
)

patterns <- c(
  'Classical Monocytes' = 'none', # Yellow-Green  
  'Intermediate Monocytes' = 'none', # Dark Orange  
  'Non-classical Monocytes' = 'none', # Medium Orchid  
  'CD15+ myeloid cells' = 'none', 
  "ARG1high CD15+" = 'stripe', # Dark Gray  
  'cDCs' = 'none')  

prop_data$new_labels <- factor(
  prop_data$new_labels,
  levels = c("ARG1high CD15+", "CD15+ myeloid cells", "Classical Monocytes", "Intermediate Monocytes", 
             "Non-classical Monocytes", "cDCs")  # von unten nach oben
)

ggplot(prop_data, aes(x = Age_Group, y = prop, fill = new_labels, pattern = new_labels)) +
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
  scale_fill_manual(values = fill_colors) +   # use colorblind-safe palette here
  scale_pattern_manual(values = patterns) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = c(0, 0)
  ) +
  labs(
    y = "Proportion (%)",
    x = "",
    fill = "Cell Type"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey80", size = 0.3),
    legend.position = "right",
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 8)
  )