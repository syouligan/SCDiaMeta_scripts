#!/bin/bash

# Load R 3.6.0
module load briglo/R/3.6.0

# Set paths
projectID='SCDiaMeta'
resultsPath='/share/ScratchGeneral/scoyou/projects/'
scriptsPath='/home/scoyou/projects/'

# Make log path
mkdir -p $scriptsPath$projectID'/'$projectID'_scripts/logs'

qsub -P OsteoporosisandTranslationalResearch -N 'SCEbuild_'$projectID -b y -wd $resultsPath$projectID'/logs' -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Make_SCE_object_from_raw_data.R'
qsub -P OsteoporosisandTranslationalResearch -N 'Prefiltering_'$projectID -b y -hold_jid 'SCEbuild_'$projectID -wd $resultsPath$projectID'/logs' -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Prefiltering_whole_experiment.R'
#qsub -P OsteoporosisandTranslationalResearch -N 'NormyCluster_'$projectID -b y -hold_jid 'Prefiltering_'$projectID -wd $resultsPath$projectID'/logs' -j y -R y -l mem_requested=8G -pe smp 64 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Whole_experiment_batch_correction_and_clustering.R'
#qsub -P OsteoporosisandTranslationalResearch -N 'Cell_Cycle_'$projectID -b y -hold_jid 'NormyCluster_'$projectID -wd $resultsPath$projectID'/logs' -j y -R y -l mem_requested=8G -pe smp 64 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Whole_experiment_cell_cycle_annotation.R'
#qsub -P OsteoporosisandTranslationalResearch -N 'WithinTissue_'$projectID -b y -hold_jid 'Cell_Cycle_'$projectID -wd $resultsPath$projectID'/logs' -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Within_tissue_batch_correction_and_clustering.R'
#qsub -P OsteoporosisandTranslationalResearch -N 'MAGIC_'$projectID -b y -hold_jid 'Cell_Cycle_'$projectID -wd $resultsPath$projectID'/logs' -j y -R y -l mem_requested=8G -pe smp 64 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/MAGIC_imputation.R'
