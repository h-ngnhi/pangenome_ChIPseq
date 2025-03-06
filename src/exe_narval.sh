#!/bin/bash 
#SBATCH --job-name=hprc_map_146
#SBATCH --output=slurm-%j.out  # %j will be replaced with the job ID
#SBATCH --error=slurm-%j.err   # %j will be replaced with the job ID
#SBATCH --account=def-bourqueg
#SBATCH --time=8:00:00
#SBATCH --mem=249G
#SBATCH --cpus-per-task=64


source pangenome_ChIPseq/src/graphs.sh
source pangenome_ChIPseq/src/genpipes.sh

# Enable xtrace to log the commands to stderr
set -x

start_time=$(date +%s)

# Set TMPDIR because the default /tmp is too small
export TMPDIR=/lustre07/scratch/hoangnhi/temp

export wd=/home/hoangnhi/projects/def-bourqueg/hoangnhi/ZNF146-507-Analysis-on-Pangenome

# Construct graphs
# Done... Will add more code if needed

# Alignment
markname=146
data_dir=$wd/Graph_genome_data
forward_trm="$data_dir/146Rep1.trim.pair1.fastq.gz $data_dir/146Rep2.trim.pair1.fastq.gz"
reverse_trm="$data_dir/146Rep1.trim.pair2.fastq.gz $data_dir/146Rep2.trim.pair2.fastq.gz"
forward_ctl="$data_dir/146Input.trim.pair1.fastq.gz"
reverse_ctl="$data_dir/146Input.trim.pair2.fastq.gz"

pipeline=vg_giraffe         # vg_giraffe or vg_map
ref=hprc-v1.1-mc-chm13
results_dir=$wd/results/${markname}_${pipeline}_${ref}
alignment $pipeline $ref $results_dir "$forward_trm" "$reverse_trm" "treatment"
alignment $pipeline $ref $results_dir "$forward_ctl" "$reverse_ctl" "control"

end_time=$(date +%s)

elapsed=$(( (end_time - start_time) / 3600 ))

# Log the total runtime
echo "Total runtime: ${elapsed} hours" >> /home/nhi/projects/stdout/${SLURM_JOB_ID}.e${SLURM_JOB_NAME}