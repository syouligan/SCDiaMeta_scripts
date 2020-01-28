#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Load Ryans preconstructed dataset, make a SingleCellExperiment Object and practice dataset
# --------------------------------------------------------------------------

# Set working directory, load data and libraries
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  location <- "/Users/mac/cloudstor/"
} else {
  location <- "/share/ScratchGeneral/scoyou/"
}

setwd(paste0(location, "projects/SCDiaMeta/Ryan_data"))

#Libraries
library('Seurat')
library('scater')
library('grid')
library('ggplot2')
library('ggridges')
library('viridis')
library('Matrix')
library('dplyr')
library('tidyr')

# Load Seurat object made by Ryan (>22Gb!)
load("Analysis_1_object_post_UMAP_res_1_260319_v3.Robj", verbose=TRUE)

# Convert Seurat object made by Ryan (>22Gb!) to SingleCellExperiement object (~1.2Gb) and save
diameta_filtered_exp <- as.SingleCellExperiment(regions.combined.final1.v3)
saveRDS(diameta_filtered_exp, "Analysis_1_object_post_UMAP_res_1_260319_v3_SCE.rds")

# Convert SingleCellExperiement object back to Seurat object (~1.4Gb) and save. No idea why it is so much smaller than Ryans Seurat object
diameta_filtered_seurat <- as.Seurat(diameta_filtered_exp)
saveRDS(diameta_filtered_seurat, "Analysis_1_object_post_UMAP_res_1_260319_v3_Seurat.rds")

# Makes practice dataset containing 5% of cells per sample. Used to develop scripts locally.
diameta_filtered_exp$cellIDs <- rownames((colData(diameta_filtered_exp)))
fe_subset <- data.frame(data.frame(colData(diameta_filtered_exp)) %>%
                          group_by(SampleID) %>%
                          sample_frac(0.05))
diameta_filtered_exp$Practice_subset <- is.element(rownames(colData(diameta_filtered_exp)), fe_subset$cellIDs)
practice_diameta_filtered_exp <- diameta_filtered_exp[,which(diameta_filtered_exp$Practice_subset)]
saveRDS(practice_diameta_filtered_exp, "Analysis_1_object_post_UMAP_res_1_260319_v3_SCE_Practice.rds")
