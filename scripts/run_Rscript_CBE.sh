#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --qos=medium
#SBATCH --partition=c
#SBATCH --mem=120G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --job-name="Rscript"
#SBATCH --output="logs/slurm-%x_%j.out"
#SBATCH --error="logs/slurm-%x_%j.err"

mkdir -p logs/

ml load build-env/f2022
ml load r/4.1.2-foss-2021b

#Rscript script_regressOut.nCount_RNA.R
#Rscript test_script_regressOut.nCount_RNA.R
#Rscript script_FindAllMarker.R
#Rscript script_DiffusionMap.R
#Rscript  script_tradeSeq.R
#Rscript script_cellcycle_v2.R # ~24h and >70G
Rscript script_cellcycle_v3.R
