#!/bin/bash 
#SBATCH --job-name=hprc_map_146
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

# Access MAPQ GenPipes
# chipseq_gp ipsc_K27 $wd/results/iPSC_K27 "-s 4-15"
# prim_bam results/iPSC_K27/linear_chm13/alignment/ipsc_K27/ipsc_K27/ipsc_K27.ipsc_K27.sorted.dup.bam
# align_stats results/K27_FLU/linear_chm13/alignment/H3K27AC_treatment1/H3K27AC/H3K27AC_treatment1.H3K27AC.sorted.dup.filtered.cleaned.bam

# bed_path=results/507/linear_chm13/peak_call/507/ZNF507  # to locate the bed file
# bed_file=507.ZNF507_peaks   # to take the file name and create downstream files
# csv_file=results/507/linear_chm13/enrichment.csv   # to write the enrichment csv file
# ref=t2t        # to use which ref genome
# blacklist=Genome/Blacklist/t2t.excluderanges.bed  # to use which blacklist file
# row_title=507_linear_chm13   # to use which row title
# enrichment_calc $bed_path $bed_file $csv_file $ref $blacklist $row_title

# Get input data
input() {
    local markname=$1
    if [ $markname == "K27_FLU" ]; then
        data_dir=$wd/results/K27_FLU/linear_chm13/trim
        forward_trm="$data_dir/H3K27AC_treatment1/H3K27AC/Treatment1.trim.pair1.fastq.gz"
        reverse_trm="$data_dir/H3K27AC_treatment1/H3K27AC/Treatment1.trim.pair2.fastq.gz"
        forward_ctl="$data_dir/H3K27AC_Input/Input/control.trim.pair1.fastq.gz"
        reverse_ctl="$data_dir/H3K27AC_Input/Input/control.trim.pair2.fastq.gz"
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
}

# Parameters
markname=iPSC_K27            # K27_FLU / 146 / 507
pipeline=vg_giraffe       # vg_giraffe or vg_map
ref=chm13                 # chm13 / L1_vcfbub / vcfbub / hprc-v1.1-mc-chm13
steps="7"         # Steps of the pipeline
mapq_edit="mapq60"  # mapq60 / ""
# Steps of the pipeline

    # 1. Construct the graph
    # 2. Split the graph
    # 3. Align the reads
    # 4. Calculate parameters for peak calling
    # 5. Split json to chromosomes
    # 6. Call peaks
    # NOTE: Haven't consider vg_map_convert hprc-v1.1-mc-chm13

chipseq_graph $markname $pipeline $ref "$steps" $mapq_edit

end_time=$(date +%s)
elapsed=$(( (end_time - start_time) / 3600 ))

# SIF_IMAGE="$wd/tools/gp.sif"
# BIND_DIR="$wd:/mnt"
# apptainer exec --contain --cleanenv --bind $BIND_DIR $SIF_IMAGE graph_peak_caller callpeaks_whole_genome_from_p_values 
# Log the total runtime
echo "Total runtime: ${elapsed} hours" 