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
    local type="$2" # linear or graph
    if [ $type == "graph" ]; then # data is from trimmed linear results
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
        elif [[ $markname == "SOX6" ]]; then
            data_dir=$wd/results/$markname/linear_chm13/trim/$markname
            forward_trm="$data_dir/$markname/treatment_1.trim.pair1.fastq.gz $data_dir/$markname/treatment_2.trim.pair1.fastq.gz"
            reverse_trm="$data_dir/$markname/treatment_1.trim.pair2.fastq.gz $data_dir/$markname/treatment_2.trim.pair2.fastq.gz"
            forward_ctl="$data_dir/Input/control_1.trim.pair1.fastq.gz $data_dir/Input/control_2.trim.pair1.fastq.gz"
            reverse_ctl="$data_dir/Input/control_1.trim.pair2.fastq.gz $data_dir/Input/control_2.trim.pair2.fastq.gz"
        elif [ $markname == "iPSC_K27" ]; then
            data_dir=$wd/data/iPSC
            forward_trm="$data_dir/p_N_iPSC_K27ac_r1_S5_R1_001.fastq.gz $data_dir/p_N_iPSC_K27ac_r2_S20_R1_001.fastq.gz"
            reverse_trm="$data_dir/p_N_iPSC_K27ac_r1_S5_R2_001.fastq.gz $data_dir/p_N_iPSC_K27ac_r2_S20_R2_001.fastq.gz"
            forward_ctl="$data_dir/p_N_iPSC_Input_r1_S18_R1_001.fastq.gz"
            reverse_ctl="$data_dir/p_N_iPSC_Input_r1_S18_R2_001.fastq.gz"
        elif [[ $markname == "s4_iPSC_K27" || $markname == "s4_naive_iPSC_K27" || $markname == "s4_iPSC_K27me3" ]]; then
            data_dir=$wd/results/$markname/linear_chm13/trim/$markname
            forward_trm="$data_dir/$markname/treatment_1.trim.pair1.fastq.gz $data_dir/$markname/treatment_2.trim.pair1.fastq.gz"
            reverse_trm="$data_dir/$markname/treatment_1.trim.pair2.fastq.gz $data_dir/$markname/treatment_2.trim.pair2.fastq.gz"
            forward_ctl="$data_dir/Input/control_1.trim.pair1.fastq.gz"
            reverse_ctl="$data_dir/Input/control_1.trim.pair2.fastq.gz"
        elif [ $markname == "YY1_2102Ep" ]; then
            data_dir=$wd/results/$markname/linear_chm13/trim/$markname
            forward_trm="$data_dir/$markname/treatment_1.trim.pair1.fastq.gz"
            reverse_trm="$data_dir/$markname/treatment_1.trim.pair2.fastq.gz"
            forward_ctl="$data_dir/Input/control_1.trim.pair1.fastq.gz"
            reverse_ctl="$data_dir/Input/control_1.trim.pair2.fastq.gz"
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
    else
        # Must edit when having time
        # if [ $markname == "K27_FLU" ]; then
        #     data_dir=$wd/results/K27_FLU/linear_chm13/trim
        #     forward_trm="$data_dir/H3K27AC_treatment1/H3K27AC/Treatment1_15.trim.pair1.fastq.gz"
        #     reverse_trm="$data_dir/H3K27AC_treatment1/H3K27AC/Treatment1_15.trim.pair2.fastq.gz"
        #     forward_ctl="$data_dir/H3K27AC_Input/Input/control_15.trim.pair1.fastq.gz"
        #     reverse_ctl="$data_dir/H3K27AC_Input/Input/control_15.trim.pair2.fastq.gz"
        # elif [[ $markname == "146" || $markname == "507" ]]; then
        #     data_dir=$wd/results/$markname/linear_chm13/trim/$markname
        #     forward_trm="$data_dir/ZNF$markname/${markname}Rep1.trim.pair1.fastq.gz $data_dir/ZNF$markname/${markname}Rep2.trim.pair1.fastq.gz"
        #     reverse_trm="$data_dir/ZNF$markname/${markname}Rep1.trim.pair2.fastq.gz $data_dir/ZNF$markname/${markname}Rep2.trim.pair2.fastq.gz"
        #     forward_ctl="$data_dir/Input/${markname}Input.trim.pair1.fastq.gz"
        #     reverse_ctl="$data_dir/Input/${markname}Input.trim.pair2.fastq.gz"
        # fi
        if [ $markname == "SOX6" ]; then
            data_dir=$wd/data/SOX6_K562
            forward_trm="$data_dir/11_ENCFF002DPJ.fastq.gz $data_dir/21_ENCFF002DUA.fastq.gz"
            reverse_trm="$data_dir/12_ENCFF002EGF.fastq.gz $data_dir/22_ENCFF002DPM.fastq.gz"
            forward_ctl="$data_dir/Input11_ENCFF002EFF.fastq.gz $data_dir/Input21_ENCFF002EFD.fastq.gz"
            reverse_ctl="$data_dir/Input12_ENCFF002EFH.fastq.gz $data_dir/Input22_ENCFF002EFA.fastq.gz"
        elif [ $markname == "s4_iPSC_K27" ]; then
            data_dir=$wd/data/new_s4_data
            forward_trm="$data_dir/p_N_iPSC_K27ac_r1_R1_001.fastq.gz $data_dir/p_N_iPSC_K27ac_r2_R1_001.fastq.gz"
            reverse_trm="$data_dir/p_N_iPSC_K27ac_r1_R2_001.fastq.gz $data_dir/p_N_iPSC_K27ac_r2_R2_001.fastq.gz"
            forward_ctl="$data_dir/p_N_iPSC_Input_r1_R1_001.fastq.gz"
            reverse_ctl="$data_dir/p_N_iPSC_Input_r1_R2_001.fastq.gz"
        elif [ $markname == "s4_iPSC_K27me3" ]; then
            data_dir=$wd/data/new_s4_data
            forward_trm="$data_dir/p_N_iPSC_K27me3_r1_R1_001.fastq.gz $data_dir/p_N_iPSC_K27me3_r2_R1_001.fastq.gz"
            reverse_trm="$data_dir/p_N_iPSC_K27me3_r1_R2_001.fastq.gz $data_dir/p_N_iPSC_K27me3_r2_R2_001.fastq.gz"
            forward_ctl="$data_dir/p_N_iPSC_Input_r1_R1_001.fastq.gz"
            reverse_ctl="$data_dir/p_N_iPSC_Input_r1_R2_001.fastq.gz"
        elif [ $markname == "s4_naive_iPSC_K27" ]; then
            data_dir=$wd/data/new_s4_data
            forward_trm="$data_dir/r_N_naive_K27ac_r1_R1_001.fastq.gz $data_dir/r_N_naive_K27ac_r2_R1_001.fastq.gz"
            reverse_trm="$data_dir/r_N_naive_K27ac_r1_R2_001.fastq.gz $data_dir/r_N_naive_K27ac_r2_R2_001.fastq.gz"
            forward_ctl="$data_dir/r_N_naive_Input_r1_R1_001.fastq.gz"
            reverse_ctl="$data_dir/r_N_naive_Input_r1_R2_001.fastq.gz"
        elif [ $markname == "YY1_2102Ep" ]; then
            data_dir=$wd/data/YY1_2102Ep
            forward_trm="$data_dir/ERR10313856_1.fastq.gz"
            reverse_trm="$data_dir/ERR10313856_2.fastq.gz"
            forward_ctl="$data_dir/ERR10313855_1.fastq.gz"
            reverse_ctl="$data_dir/ERR10313855_2.fastq.gz"
        fi
    fi
}

# Parameters
# Define the parameters
mark="$1"
pipeline="$2"
ref="$3"
steps="$4"
mapq_troubleshoot="$5"
giraffe_extraoptions="$6"
mapq_filter="$7"
k="$8"
w="$9"
# echo "Running chipseq_graph for mark: $this_mark"
chipseq_graph "$mark" "$pipeline" "$ref" "$steps" "$mapq_troubleshoot" "$giraffe_extraoptions" "$mapq_filter" "$k" "$w"

##########################
mark="YY1_2102Ep"
ref="" # empty means chm13 from mugqic
steps="3" 
    # Steps of the pipeline
            # 1. Create readset, design files
            # 2. Create ini file
            # 3. Run GenPipes
extra_options="" # extra options for GenPipes
# chipseq_linear "$mark" "$ref" "$steps" "$extra_options"

end_time=$(date +%s)
elapsed=$(( (end_time - start_time) / 3600 ))

# Log the total runtime
echo "Total runtime: ${elapsed} hours" 
