#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Collate all Diaphasis, Metaphasis and Marrow data into SingleCellExperiment object.
# --------------------------------------------------------------------------

# Set working directory, load data and libraries
# --------------------------------------------------------------------------

# Working directory
# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  location <- "/Users/mac/cloudstor/"
  place <- "local"
} else {
  location <- "/share/ScratchGeneral/scoyou/"
  place <- "wolfpack"
}

setwd(paste0(location, "projects/SCDiaMeta/scott_data"))

# Libraries
library('DropletUtils')
library('dplyr')
library('tidyr')
library('scater')
library('Matrix')
library('SingleCellExperiment')
library('Seurat')

# Load all data into single object
# --------------------------------------------------------------------------

# Load all samples (less samples where pooling IDs were convoluted)
data_directory <- "/share/ClusterShare/thingamajigs/wenkho/projects/Ryan/NovaSeq/181024"
raw_experiment <- read10xCounts(paste0(list.files(data_directory, full.names = TRUE)), col.names = TRUE)
raw_experiment$Sample <- gsub("/share/ClusterShare/thingamajigs/wenkho/projects/Ryan/NovaSeq/181024/", "",raw_experiment$Sample)

# Load two samples that were deconvoluted and convert to SCEs.
load("/share/ScratchGeneral/ryacha/projects/Wellcome/single_cell/Diaphysis_Metaphysis_Marrow_merged_minus_marrow_2_newly_M1_D4_deconvoluted/metaphysis_1_and_diaphysis_4_decon_190319/metaphysis_1_decon_read10x_object_190319.Robj", verbose = TRUE)
load("/share/ScratchGeneral/ryacha/projects/Wellcome/single_cell/Diaphysis_Metaphysis_Marrow_merged_minus_marrow_2_newly_M1_D4_deconvoluted/metaphysis_1_and_diaphysis_4_decon_190319/diaphysis_4_decon_read10x_object_190319.Robj", verbose = TRUE)

meta1_sce <- SingleCellExperiment(assays = list(counts = Metaphysis_1))
meta1_sce$Sample <- "Metaphysis_1"
if(identical(gsub('\\..*', "", rownames(meta1_sce)), gsub('\\..*', "", rowData(raw_experiment)$Symbol))) {
rownames(meta1_sce) <- rownames(raw_experiment)
rowData(meta1_sce)$ID <- rowData(raw_experiment)$ID
rowData(meta1_sce)$Symbol <- rowData(raw_experiment)$Symbol
} else {
  print("Not identical. Check whats going on.")
  }
meta1_sce$Barcode <- paste0(colnames(meta1_sce), "-1")
colnames(meta1_sce) <- paste0("14_", colnames(meta1_sce), "-1")

dia4_sce <- SingleCellExperiment(assays = list(counts = Diaphysis_4))
dia4_sce$Sample <- "Diaphysis_4"
if(identical(gsub('\\..*', "", rownames(dia4_sce)), gsub('\\..*', "", rowData(raw_experiment)$Symbol))) {
  rownames(dia4_sce) <- rownames(raw_experiment)
  rowData(dia4_sce)$ID <- rowData(raw_experiment)$ID
  rowData(dia4_sce)$Symbol <- rowData(raw_experiment)$Symbol
} else {
  print("Not identical. Check whats going on.")
}
dia4_sce$Barcode <- paste0(colnames(dia4_sce), "-1")
colnames(dia4_sce) <- paste0("15_", colnames(dia4_sce), "-1")

# Merge all samples into a single SCE object
combined.sce <- BiocGenerics::cbind(raw_experiment, meta1_sce, dia4_sce, deparse.level=1)

# Annotate cells with sample information
sample_IDs <- unique(combined.sce$Sample)
sample_IDs <- data.frame(sample_IDs) %>%
  separate(sample_IDs, c("Tissue", "Replicate"), "_", remove = FALSE)
idx <- match(combined.sce$Sample, sample_IDs$sample_IDs)
combined.sce$Tissue <- as.character(sample_IDs$Tissue) [idx]
combined.sce$Replicate <- as.character(sample_IDs$Replicate) [idx]

# Make unique gene names by combining Ensembl ID and geneSymbol
rowData(combined.sce)$Ensembl <- as.character(rowData(combined.sce)$ID)
rowData(combined.sce)$GeneSymbol <- as.character(gsub('\\..*', "", rowData(combined.sce)$Symbol))
rownames(combined.sce) <- uniquifyFeatureNames(rowData(combined.sce)$Ensembl, rowData(combined.sce)$GeneSymbol)

# Save experiment before filtering
saveRDS(combined.sce, "DiaMeta_dataset_prefiltering.rds")
