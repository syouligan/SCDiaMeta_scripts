#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Collate all Diaphasis, Metaphasis and Marrow data into SingleCellExperiment object.
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

setwd(paste0(location, "projects/SCDiaMeta/project_results/prefiltered"))

# Libraries
library('DropletUtils')
library('dplyr')
library('tidyr')
library('scater')
library('EnsDb.Mmusculus.v79')
library('grid')
library('ggplot2')
library('ggridges')
library('viridis')
library('Matrix')

set.seed(100)
# Load unfiltered SCE object
combined.sce <- readRDS(paste0(location, "projects/SCDiaMeta/scott_data/DiaMeta_dataset_prefiltering.rds"))

# Identifiy and remove empty Droplets. NOTE: you dont need to do this for Cell Ranger pipelines v3 or later (v2 has a genes.tsv, v3 has a features.tsv file in Cell ranger results folder)
empty.out <- emptyDrops(counts(combined.sce))
summary(empty.out$FDR <= 0.001)
table(Sig = empty.out$FDR <= 0.001, Limited = empty.out$Limited)
filtered_exp <- combined.sce[,which(empty.out$FDR <= 0.001)]

# Idenitify cells to discard based on 3MAD outlier in either number of detected genes, library size or mitochondrial content
location <- mapIds(EnsDb.Mmusculus.v79, keys=rowData(filtered_exp)$Ensembl, column="SEQNAME", keytype="GENEID")
stats <- perCellQCMetrics(filtered_exp, subsets=list(Mito=which(location=="MT")))
filtered_exp$Lib_size <- stats$sum
filtered_exp$Genes_detected <- stats$detected
filtered_exp$Mito_percent <- stats$subsets_Mito_percent

discard <- quickPerCellQC(stats, percent_subsets=c("subsets_Mito_percent"), batch=filtered_exp$Sample)
filtered_exp$discard_LibGenes <- discard$low_lib_size | discard$low_n_features
filtered_exp$discard_Mito <- discard$high_subsets_Mito_percent
filtered_exp$discard <- discard$discard
discard_stats <- DataFrame(colSums(as.matrix(discard)))
colnames(discard_stats) <- "Cell#"
print(discard_stats)

# Removes cells based on hard filters 
# filtered_exp$discard_Genes <- discard$Genes_detected < 300
# filtered_exp$discard_Mito <- discard$Mito_percent > 20
# filtered_exp$discard_Lib <- discard$Lib_size < 1000
# filtered_exp$discard <- filtered_exp$discard_Genes | filtered_exp$discard_Mito | filtered_exp$discard_Lib
# 
# sum(discard$discard_Genes)
# sum(discard$discard_Mito)
# sum(discard$discard_Lib)

# Plot QC stats
ggplot(data.frame(colData(filtered_exp)), aes(x = Lib_size, y = Sample, fill = Tissue)) +
  geom_density_ridges() +
  theme_minimal() +
  ggsave("Library_size_ridge_before.pdf", useDingbats = FALSE)

ggplot(data.frame(colData(filtered_exp)), aes(x = Genes_detected, y = Sample, fill = Tissue)) +
  geom_density_ridges() +
  theme_minimal() +
  ggsave("Number_of_genes_ridge_before.pdf", useDingbats = FALSE)

ggplot(data.frame(colData(filtered_exp)), aes(x = Mito_percent, y = Sample, fill = Tissue)) +
  geom_density_ridges() +
  theme_minimal() +
  ggsave("Mito_percent_ridge_before.pdf", useDingbats = FALSE)

# Remove "discard" cells
filtered_exp <- filtered_exp[ ,which(!filtered_exp$discard)]

ggplot(data.frame(colData(filtered_exp)), aes(x = Lib_size, y = Sample, fill = Tissue)) +
  geom_density_ridges() +
  theme_minimal() +
  ggsave("Library_size_ridge_after.pdf", useDingbats = FALSE)

ggplot(data.frame(colData(filtered_exp)), aes(x = Genes_detected, y = Sample, fill = Tissue)) +
  geom_density_ridges() +
  theme_minimal() +
  ggsave("Number_of_genes_ridge_after.pdf", useDingbats = FALSE)

ggplot(data.frame(colData(filtered_exp)), aes(x = Mito_percent, y = Sample, fill = Tissue)) +
  geom_density_ridges() +
  theme_minimal() +
  ggsave("Mito_percent_ridge_after.pdf", useDingbats = FALSE)

# Remove genes without counts in at least 3 cells in each samples for any tissue
tmp_structure <- data.frame(unique(colData(filtered_exp)[ ,c("Tissue", "Sample")]))
tmp_table <- data.frame(table(tmp_structure$Tissue))
colnames(tmp_table) <- c("Tissue", "Reps")

GOI <- c()
for(t in tmp_table$Tissue) {
  tissue_exp <- filtered_exp[,filtered_exp$Tissue == t]
  sampleGT3 <- c()
  for(i in unique(tissue_exp$Sample)) {
    tmp_sample <- tissue_exp[,tissue_exp$Sample == i]
    sampleGT3 <- cbind(sampleGT3, rowSums(counts(tmp_sample) > 0) > 3) # Greater than 0 counts in at least 3 cells in sample
  }
  GOI <- cbind(GOI, rowSums(sampleGT3) == tmp_table[tmp_table$Tissue == t, "Reps", drop = TRUE]) # Above 3 in each replicate of each tissue type
}
colnames(GOI) <- paste0(tmp_table$Tissue, "_active")
GOI <- data.frame(GOI)
GOI$Any_Active <- rowSums(GOI) > 0 # Above three in any replicate of any tissue type
rowData(filtered_exp) <- cbind(rowData(filtered_exp), GOI) # Add to rowData
filtered_exp <- filtered_exp[which(GOI$Any_Active), ] # Subset to active genes

colSums(GOI) # Number of genes "active" in each tissue

# Number of cells remaining per sample
cells_remaining <- data.frame(table(filtered_exp$Sample))
write.csv(cells_remaining, "Cells_remaining_number.csv", row.names = FALSE)

ggplot(data=cells_remaining, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", color="black") +
  ggtitle("Cells remaining") +
  scale_color_viridis(discrete=TRUE) +
  coord_flip() +
  theme_classic() +
  ggsave("Cells_remaining_barplot.pdf", useDingbats = FALSE)

# Save datasets.
# --------------------------------------------------------------------------

# Save practice dataset (5% of cells from each sample)
filtered_exp$cellIDs <- rownames((colData(filtered_exp)))
fe_subset <- data.frame(data.frame(colData(filtered_exp)) %>%
                          group_by(Sample) %>%
                          sample_frac(0.05))
filtered_exp$Practice_subset <- is.element(rownames(colData(filtered_exp)), fe_subset$cellIDs)
practice_exp <- filtered_exp[,which(filtered_exp$Practice_subset)]

if (place == "wolfpack") {
  saveRDS(practice_exp, "practice_all_data/Prefiltered_experiment_Practice.rds")
} else {
  print("Working local")
}

combined.sce$cellIDs <- rownames((colData(combined.sce)))
fe_subset <- data.frame(data.frame(colData(combined.sce)) %>%
                          group_by(Sample) %>%
                          sample_frac(0.05))
combined.sce$Practice_subset <- is.element(rownames(colData(combined.sce)), fe_subset$cellIDs)
practice_exp <- combined.sce[,which(combined.sce$Practice_subset)]
if (place == "wolfpack") {
  saveRDS(practice_exp, "practice_all_data/combined.sce_Practice.rds")
} else {
  print("Working local")
}

# Save total filtered dataset
if (place == "wolfpack") {
  saveRDS(filtered_exp, "all_data/Prefiltered_experiment_All.rds")
} else {
  print("Working local")
}

# Save individual samples
if (place == "wolfpack") {
  for(i in unique(colData(filtered_exp)$Sample)) {
    dir.create(paste0("individual/", i))
    saveRDS(filtered_exp[,filtered_exp$Sample == i], paste0("individual/", i,"/Prefiltered_experiment_", i, ".rds"))
  }
} else {
  print("Working local")
}