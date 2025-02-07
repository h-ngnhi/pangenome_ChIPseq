#!/bin/bash 
source pangenome_ChIPseq/src/graph_chipseq_ana.sh

#SBATCH --job-name=vg_giraffe_construct pangenome_ChIPseq/src/exe_narval.sh
#SBATCH --output=slurm-%j.out  # %j will be replaced with the job ID
#SBATCH --error=slurm-%j.err   # %j will be replaced with the job ID
#SBATCH --time=8:00:00
#SBATCH --mem=249G
#SBATCH --cpus-per-task=64
