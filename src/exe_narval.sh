#!/bin/bash 
#SBATCH --job-name=hprc_map_146
#SBATCH --output=slurm-%j.out  # %j will be replaced with the job ID
#SBATCH --error=slurm-%j.err   # %j will be replaced with the job ID
#SBATCH --account=def-bourqueg
#SBATCH --time=24:00:00
#SBATCH --mem=249G
#SBATCH --cpus-per-task=64


source pangenome_ChIPseq/src/graphs.sh
source pangenome_ChIPseq/src/genpipes.sh

# Enable xtrace to log the commands to stderr
set -x

start_time=$(date +%s)

# Set TMPDIR because the default /tmp is too small
export TMPDIR=/lustre07/scratch/hoangnhi/temp

export wd=$(pwd)

# Construct graphs
# Done... Will add more code if needed
# vg_map_convert hprc-v1.1-mc-chm13

markname=146
pipeline=vg_giraffe         # vg_giraffe or vg_map
ref=chm13                    # chm13 / L1_vcfbub / vcfbub / hprc-v1.1-mc-chm13
results_dir=$wd/results/${markname}/${pipeline}_${ref}

input() {
    local markname=$1
    if [ $markname == "K27_FLU" ]; then
        data_dir=$wd/cgroza_data/H3K27AC_FLU
        forward_trm="$data_dir/treatment/H3K27AC.forward_treatment_1.fastq.gz"
        reverse_trm="$data_dir/treatment/H3K27AC.reverse_treatment_1.fastq.gz"
        forward_ctl="$data_dir/control/H3K27AC.forward_control.fastq.gz"
        reverse_ctl="$data_dir/control/H3K27AC.reverse_control.fastq.gz"
        trm_json="treatment_alignments.filtered.json"
        ctl_json=""
    elif [[ $markname == "146" || $markname == "507" ]]; then
        data_dir=$wd/results/$markname/linear_chm13/trim/$markname/ZNF$markname
        forward_trm="$data_dir/${markname}Rep1.trim.pair1.fastq.gz $data_dir/${markname}Rep2.trim.pair1.fastq.gz"
        reverse_trm="$data_dir/${markname}Rep1.trim.pair2.fastq.gz $data_dir/${markname}Rep2.trim.pair2.fastq.gz"
        forward_ctl="$data_dir/${markname}Input.trim.pair1.fastq.gz"
        reverse_ctl="$data_dir/${markname}Input.trim.pair2.fastq.gz"
        trm_json="treatment_alignments.filtered.json"
        ctl_json="control_alignments.filtered.json"
    fi
}

input $markname

# Alignment
# alignment $pipeline $ref $results_dir "$forward_trm" "$reverse_trm" "treatment"
# if [ -n "$forward_ctl" ]; then
#     alignment $pipeline $ref $results_dir "$forward_ctl" "$reverse_ctl" "control"
# fi

#################
# Peak calling
callpeaks "$trm_json" "$ctl_json" "${results_dir#$wd/}" "$ref" $pipeline "narval" $wd/tools/gp.sif


end_time=$(date +%s)

elapsed=$(( (end_time - start_time) / 3600 ))

# Log the total runtime
echo "Total runtime: ${elapsed} hours" 