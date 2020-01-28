#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Examine CCR expression across cell populations in the bone marrow.
# --------------------------------------------------------------------------

# Set working directory, load data and libraries
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  location <- "/Users/mac/cloudstor/"
  place <- "local"
} else {
  location <- "/share/ScratchGeneral/scoyou/"
  place <- "wolfpack"
}

setwd(paste0(location, "projects/SCDiaMeta/project_results/chemokines"))

# Libraries
library('DropletUtils')
library('dplyr')
library('tidyr')
library('scater')
library('EnsDb.Hsapiens.v75')
library('grid')
library('ggplot2')
library('ggridges')
library('viridis')
library('Matrix')
library('reshape2')

# Read CSV with Chemokine Ligand-Receptor detail
CCRs <- read.csv("20200124_Chemokine_receptor_targets_curated.csv", header = TRUE)

# Read CSV with cel type annotation for clusters
cell_types <- read.csv(paste0(location, "projects/SCDiaMeta/Ryan_data/Cell_type_classification_whole_bone_res_1_210619.csv"), header = TRUE)

# Make single cell object using counts matrix
if (place == "local") {
  
  filtered_exp <- readRDS(paste0(location, "projects/SCDiaMeta/Ryan_data/Analysis_1_object_post_UMAP_res_1_260319_v3_SCE_Practice.rds"))
  
} else {

  filtered_exp <- readRDS(paste0(location, "projects/SCDiaMeta/Ryan_data/Analysis_1_object_post_UMAP_res_1_260319_v3_SCE.rds"))

  }

# Add cell type IDs to SCE object
filtered_exp$Ryan_cluster <- as.character(filtered_exp$ident)
cell_type_key <- as.character(cell_types$Cell_type)
names(cell_type_key) <- as.character(cell_types$Cluster_ID)
cell_type_vec <- filtered_exp$Ryan_cluster
for(i in as.character(unique(cell_type_vec))){
  cell_type_vec[cell_type_vec == i] <- cell_type_key[[i]]
  } # Makes vector of colors according to tissue
filtered_exp$Cell_types <- cell_type_vec

# Subset to single cells and chemokine receptor genes
CCR.sce <- filtered_exp[unique(as.character(CCRs$Mouse_receptor_GeneSymbol)), ]

# For each cell type calculate the number of cells positive for CCRs
CCR.percent <- c()
for (i in unique(filtered_exp$Cell_types)) {
  CCR.cell <- rowSums(counts(CCR.sce[,CCR.sce$Cell_types == i]) > 0)/ncol(CCR.sce[,CCR.sce$Cell_types == i])*100 # Calculates percent of cell types that express each receptor
  CCR.percent <- cbind(CCR.percent, CCR.cell)
}

CCR.percent <- data.frame(CCR.percent)
colnames(CCR.percent) <- unique(filtered_exp$Cell_types)
CCR.percent$Gene <- rownames(CCR.percent)

CCR.percent <- CCR.percent[order(CCR.percent$Gene),]

# Make dotplot relecting percentage of expressing cells per CCR
CCR.percent.melt <- melt(CCR.percent)
colnames(CCR.percent.melt) <- c("Gene", "Cell_type", "Percent")
CCR.percent.melt$Gene <- factor(CCR.percent.melt$Gene, levels = rev(CCR.percent$Gene))
CCR.percent.melt$Percent[CCR.percent.melt$Percent == 0] <- NA # Removes 0 values from plot if desired

ggplot(CCR.percent.melt, aes(x = Gene, y = Cell_type)) +
  geom_point(aes(size = Percent, color = Percent)) +
  scale_color_viridis_c(option = "B", limits = c(0, 100)) +
  scale_radius(range = c(0, 6), limits = c(0, 100)) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggsave("DiaMeta_CCR_expression.pdf", useDingbats = FALSE)



