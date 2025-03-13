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

markname=507
pipeline=vg_giraffe         # vg_giraffe or vg_map
ref=chm13                    # chm13 / L1_vcfbub / vcfbub / hprc-v1.1-mc-chm13
results_dir=$wd/results/${markname}/${pipeline}_${ref}

input() {
    local markname=$1
    if [ $markname == "K27_FLU" ]; then
        data_dir=$wd/cgroza_data/H3K27AC_FLU
        forward_trm="$data_dir/treatment/H3K27AC.forward_treatment_1.fastq.gz"
        reverse_trm="$data_dir/treatment/H3K27AC.reverse_treatment_1.fastq.gz"
        forward_ctl=""
        reverse_ctl=""
    elif [[ $markname == "146" || $markname == "507" ]]; then
        data_dir=$wd/results/$markname/${markname}_linear_t2t/trim/$markname/ZNF$markname
        forward_trm="$data_dir/${markname}Rep1.trim.pair1.fastq.gz $data_dir/${markname}Rep2.trim.pair1.fastq.gz"
        reverse_trm="$data_dir/${markname}Rep1.trim.pair2.fastq.gz $data_dir/${markname}Rep2.trim.pair2.fastq.gz"
        forward_ctl="$data_dir/${markname}Input.trim.pair1.fastq.gz"
        reverse_ctl="$data_dir/${markname}Input.trim.pair2.fastq.gz"
    fi
}

input $markname

# Alignment
alignment $pipeline $ref $results_dir "$forward_trm" "$reverse_trm" "treatment"
if [ -n "$forward_ctl" ]; then
    alignment $pipeline $ref $results_dir "$forward_ctl" "$reverse_ctl" "control"
fi

#################
# Peak calling
fragment_length() {
    grep -i "predicted fragment length is" | head -n1 | sed -E 's/.*predicted fragment length is[[:space:]]*([0-9]+).*/\1/'
}
module load apptainer
SIF_IMAGE="$wd/tools/macs2_latest.sif"
BIND_DIR="$results_dir:/mnt"
frag_len=$(apptainer exec --contain --cleanenv --bind $BIND_DIR $SIF_IMAGE macs2 predictd -i /mnt/treatment_alignments.bam 2>&1 | fragment_length)
read_len=$(zcat $forward_trm | head -2 | tail -1 | wc -c)
unique_reads=$(grep -Po '"sequence": "\K([ACGTNacgtn]{20,})"' $results_dir/treatment_alignments.filtered.json | sort | uniq | wc -l)
gp=$wd/tools/gp.sif
callpeaks "treatment_alignments.filtered.json" "control_alignments.filtered.json" "${results_dir#$wd/}" "$ref" $pipeline $frag_len $read_len $unique_reads $gp


end_time=$(date +%s)

elapsed=$(( (end_time - start_time) / 3600 ))

# Log the total runtime
echo "Total runtime: ${elapsed} hours" 