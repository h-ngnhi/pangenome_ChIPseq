#!/bin/bash -l 
#SBATCH --job-name=hprc_map
#SBATCH --mem=249G
#SBATCH --time=24:00:00
#SBATCH -n 1
#SBATCH -c 64
#SBATCH -o /home/nhi/projects/stdout/%j.o%x
#SBATCH -e /home/nhi/projects/stdout/%j.e%x
#SBATCH -p node
source pangenome_ChIPseq/src/graphs.sh
source pangenome_ChIPseq/src/genpipes.sh

# Enable xtrace to log the commands to stderr
set -x

start_time=$(date +%s)

# Set TMPDIR because the default /tmp is too small
# export TMPDIR=/lustre07/scratch/hoangnhi/temp #make this input for the pipeline
export TMPDIR=/home/nhi/tmp

export wd=$(pwd)

#######
#1. Construct graphs
variant_vcf=$wd/genome_data/L1_annotation.vcfbub.vcf.gz
ref=$wd/genome_data/chm13v2.0.fa
## Construct + index giraffe graph
# vg_giraffe_graph "" $ref chm13
# vg_giraffe_graph $variant_vcf $ref L1_vcfbub
## Construct + index vgmap graph
# bgzip Graph_genome_data/L1_annotation/L1_annotation.vcfbub.vcf
# tabix -p vcf Graph_genome_data/L1_annotation/L1_annotation.vcfbub.vcf.gz
# vg_map_graph "" $ref chm13
# vg_map_graph $variant_vcf $ref L1_vcfbub

#######
#2. Alignment
markname=ipsc_K27
data_dir=/home/share/saitoulab/liu_project/kumadata
forward_trm="$data_dir/p_N_iPSC_K27ac_r1_S5_R1_001.fastq.gz $data_dir/p_N_iPSC_K27ac_r2_S20_R1_001.fastq.gz"
reverse_trm="$data_dir/p_N_iPSC_K27ac_r1_S5_R1_001.fastq.gz $data_dir/p_N_iPSC_K27ac_r2_S20_R2_001.fastq.gz"
forward_ctl="$data_dir/p_N_iPSC_Input_r1_S18_R1_001.fastq.gz"
reverse_ctl="$data_dir/p_N_iPSC_Input_r1_S18_R2_001.fastq.gz"

# GenPipes
# result_dir=$wd/results/${markname}_linear_chm13
# mkdir -p $result_dir
# create_readset $markname "$forward_trm" "$reverse_trm" "$forward_ctl" "$reverse_ctl" $result_dir
# create_design $markname $result_dir
# create_config $wd/genome_data/GenPipes_t2t/Homo_sapiens.T2T-CHM13v2.0.maskedY.rCRS.EBV.fa
# chipseq_gp $markname $result_dir

# Graphs
# module load samtools
# cd results/ipsc_K27_vg_map_L1_vcfbub_converted
# samtools view treatment_alignments.bam | awk '{print $5}' > mapq_scores.txt
# awk '$1 ==0 {count++} END {print "MAPQ = 0:", count+0}' mapq_scores.txt
# awk '$1 >0 && $1 < 30 {count++} END {print "0 < MAPQ < 30:", count+0}' mapq_scores.txt
# awk '$1 >= 30 && $1 < 60 {count++} END {print "30 <= MAPQ < 60:", count+0}' mapq_scores.txt
# awk '$1 == 60 {count++} END {print "MAPQ = 60:", count+0}' mapq_scores.txt
pipeline=vg_map
ref=L1_vcfbub_converted
results_dir=$wd/results/${markname}_${pipeline}_${ref}

vg_map_convert L1_vcfbub
alignment $pipeline $ref $results_dir "$forward_trm" "$reverse_trm" "treatment"

alignment $pipeline $ref $results_dir "$forward_ctl" "$reverse_ctl" "control"

# split_graph $pipeline $ref $wd/tools/gp.sif


#################
# Report some parameter before calling peak
# # fragment length
# module load macs2
# macs2 predictd -i $wd/results/ipsc_K27_vg_map_chm13/treatment_alignments.bam
# # read length
# fastq=$data_dir/p_N_iPSC_K27ac_r1_S5_R1_001.fastq.gz
# echo $(zcat $fastq | head -2 | tail -1 | wc -c)
# # unique reads
# json=$wd/results/ipsc_K27_vg_giraffe_chm13/treatment_alignments.filtered.json
# echo $(grep -Po '"sequence": "\K([ACGTNacgtn]{20,})"' $json | sort | uniq | wc -l)

# graph_type=hprc-v1.1-mc-chm13
# pipeline=vg_giraffe
# fragment_length=120
# read_length=102
# unique_reads=13295105
# callpeaks "treatment_alignments.filtered.json" "" "${results_dir#$wd/}" "$graph_type" $pipeline $fragment_length $read_length $unique_reads


# # vg gamcompare  -T $wd/Graph_genome_chm13_trimmomatic/146_alignments.gam $wd/Graph_genome_vcfbub_trimmomatic/146_alignments.gam -t 64 > $wd/Graph_genome_data/146_alignments_compare.tsv
# awk -F'\t' '$1 == 1 {correct_variant++} $1 == 0 {incorrect_variant++} END {
#   print "Reads aligning better in the variant graph: " correct_variant;
#   print "Reads aligning better in the linear graph: " incorrect_variant;
# }' results/compare.tsv
# vg gamcompare -T results/ipsc_K27_vg_giraffe_L1_vcfbub/treatment_alignments.gam results/ipsc_K27_vg_map_L1_vcfbub_converted/treatment_alignments.gam -t 64 > results/compare.tsv

#-------------------
# vg inject to inject the GenPipe alignments into the graph then call peaks
# graph_type=L1_vcfbub
# data_dir=$wd/cgroza_data/H3K27AC_FLU
#
# mkdir -p $results_dir
# # samtools view -h cgroza_data/H3K27AC_FLU/H3K27AC_CHM13linear/alignment/H3K27AC_treatment1/H3K27AC/H3K27AC_treatment1.H3K27AC.sorted.bam | awk '$3 != "chrEBV" && $3 != "chrY" || $1 ~ /^@/' | samtools view -b -o $results_dir/genpipes_without_chrYEBV.bam
# # vg inject -x $vg_index.gbz $results_dir/genpipes_without_chrY.bam > $results_dir/treatment_alignments_genpipes.gam
# # vg view -aj $results_dir/treatment_alignments_genpipes.gam > $results_dir/treatment_alignments_genpipes.json
fragment_length=212
read_length=102
unique_reads=112938149
# split_graph $pipeline $graph_type
# callpeaks "treatment_alignments.filtered.json" "control_alignments.filtered.json" "${results_dir#$wd/}" "$ref" $pipeline $fragment_length $read_length $unique_reads

end_time=$(date +%s)

elapsed=$(( (end_time - start_time) / 3600 ))

# Log the total runtime
echo "Total runtime: ${elapsed} hours" >> /home/nhi/projects/stdout/${SLURM_JOB_ID}.e${SLURM_JOB_NAME}