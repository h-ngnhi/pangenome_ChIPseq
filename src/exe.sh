#!/bin/bash 
#SBATCH --output=slurm-%j.out  # %j will be replaced with the job ID
#SBATCH --error=slurm-%j.err   # %j will be replaced with the job ID
#SBATCH --account=def-bourqueg
#SBATCH --time=24:00:00
#SBATCH --mem=249G
#SBATCH --cpus-per-task=64


source pangenome_ChIPseq/src/graphs.sh
source pangenome_ChIPseq/src/linear.sh

# Enable xtrace to log the commands to stderr
set -x

start_time=$(date +%s)

# Set TMPDIR because the default /tmp is too small
export TMPDIR=/lustre07/scratch/hoangnhi/temp

export wd=$(pwd)

# Get input data
input() {
    local markname=$1
    # if [ -z "$k" ]; then 
        if [ $markname == "K27_FLU" ]; then
            data_dir=$wd/results/K27_FLU/linear_chm13/trim
            forward_trm="$data_dir/H3K27AC_treatment1/H3K27AC/Treatment1_15.trim.pair1.fastq.gz"
            reverse_trm="$data_dir/H3K27AC_treatment1/H3K27AC/Treatment1_15.trim.pair2.fastq.gz"
            forward_ctl="$data_dir/H3K27AC_Input/Input/control_15.trim.pair1.fastq.gz"
            reverse_ctl="$data_dir/H3K27AC_Input/Input/control_15.trim.pair2.fastq.gz"
        elif [[ $markname == "146" || $markname == "507" ]]; then
            data_dir=$wd/results/$markname/linear_chm13/trim/$markname
            forward_trm="$data_dir/ZNF$markname/${markname}Rep1.trim.pair1.fastq.gz $data_dir/ZNF$markname/${markname}Rep2.trim.pair1.fastq.gz"
            reverse_trm="$data_dir/ZNF$markname/${markname}Rep1.trim.pair2.fastq.gz $data_dir/ZNF$markname/${markname}Rep2.trim.pair2.fastq.gz"
            forward_ctl="$data_dir/Input/${markname}Input.trim.pair1.fastq.gz"
            reverse_ctl="$data_dir/Input/${markname}Input.trim.pair2.fastq.gz"
        elif [ $markname == "iPSC_K27" ]; then
            data_dir=$wd/data/iPSC
            forward_trm="$data_dir/p_N_iPSC_K27ac_r1_S5_R1_001.fastq.gz $data_dir/p_N_iPSC_K27ac_r2_S20_R1_001.fastq.gz"
            reverse_trm="$data_dir/p_N_iPSC_K27ac_r1_S5_R2_001.fastq.gz $data_dir/p_N_iPSC_K27ac_r2_S20_R2_001.fastq.gz"
            forward_ctl="$data_dir/p_N_iPSC_Input_r1_S18_R1_001.fastq.gz"
            reverse_ctl="$data_dir/p_N_iPSC_Input_r1_S18_R2_001.fastq.gz"
        fi
        # elif [ $markname == "iPSC_K27" ]; then
        #     data_dir=$wd/data/iPSC
        #     forward_trm="$data_dir/Treatment1.trim.pair1.fastq.gz $data_dir/Treatment2.trim.pair1.fastq.gz"
        #     reverse_trm="$data_dir/Treatment1.trim.pair2.fastq.gz $data_dir/Treatment2.trim.pair2.fastq.gz"
        #     forward_ctl="$data_dir/Input.trim.pair1.fastq.gz"
        #     reverse_ctl="$data_dir/Input.trim.pair2.fastq.gz"
        # fi
    # else # Use smaller dataset to test for k and w, no need to care about alignment of control
    #     if [ $markname == "iPSC_K27" ]; then
    #         data_dir=$wd/data/iPSC
    #         forward_trm="$data_dir/Treatment.pair1.chr6.42M_45M.8.fastq.gz"
    #         reverse_trm="$data_dir/Treatment.pair2.chr6.42M_45M.8.fastq.gz"
    #         forward_ctl=""
    #         reverse_ctl=""
    #     elif [ $markname == "K27_FLU" ]; then
    #         data_dir=$wd/results/K27_FLU/linear_chm13/trim
    #         forward_trm="$data_dir/H3K27AC_treatment1/H3K27AC/Treatment.pair1.chr6.42M_45M.15.fastq.gz"
    #         reverse_trm="$data_dir/H3K27AC_treatment1/H3K27AC/Treatment.pair2.chr6.42M_45M.15.fastq.gz"
    #         forward_ctl=""
    #         reverse_ctl=""
    #     fi
    # fi
}

# Parameters
# Define the parameters
mark="$1"
pipeline="$2"
ref="$3"
steps="$4"
mapq_troubleshoot="$5"
mapq_filter="$6"
k="$7"
w="$8"
# echo "Running chipseq_graph for mark: $this_mark"
chipseq_graph "$mark" "$pipeline" "$ref" "$steps" "$mapq_troubleshoot" "$mapq_filter" "$k" "$w"

end_time=$(date +%s)
elapsed=$(( (end_time - start_time) / 3600 ))

# Log the total runtime
echo "Total runtime: ${elapsed} hours" 
