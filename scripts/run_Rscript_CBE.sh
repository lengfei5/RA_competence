#!/bin/bash

#SBATCH --time=2-00:00:00
#SBATCH --qos=medium
#SBATCH --partition=c
#SBATCH --mem=120G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name="Rscript"
#SBATCH --output="logs/slurm-%x_%j.out"
#SBATCH --error="logs/slurm-%x_%j.err"

mkdir -p logs/

ml load build-env/f2022
ml load r/4.1.2-foss-2021b
ml load r-bundle-bioconductor/3.14-foss-2021b-r-4.1.2

#Rscript script_regressOut.nCount_RNA.R
#Rscript test_script_regressOut.nCount_RNA.R
#Rscript script_FindAllMarker.R
#Rscript script_DiffusionMap.R
#Rscript  script_tradeSeq.R
#Rscript script_cellcycle_v2.R # ~24h and >70G
#Rscript script_cellcycle_v3.R
#Rscript impuate_scRNA_SAVER.R
#Rscript run_multiome_merging.R
Rscript script_regressOut.nCount_multiome.R
