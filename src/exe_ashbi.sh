#!/bin/bash -l 
#SBATCH --job-name=vg_map_vcfbub_construct
#SBATCH --mem=249G
#SBATCH --time=8:00:00
#SBATCH -n 1
#SBATCH -c 64
#SBATCH -o /home/nhi/projects/stdout/%x.o%j
#SBATCH -e /home/nhi/projects/stdout/%x.e%j
#SBATCH -p node
source pangenome_ChIPseq/src/graph_chipseq_ana.sh

# Enable xtrace to log the commands to stderr
set -x

start_time=$(date +%s)

# Set TMPDIR because the default /tmp is too small
# export TMPDIR=/lustre07/scratch/hoangnhi/temp #make this input for the pipeline
export TMPDIR=/tmp

export wd=$(pwd)

#######
variant_vcf=$wd/genome_data/L1_annotation/L1_annotation.vcfbub.vcf.gz
ref=$wd/genome_data/chm13v2.0.fa
# Construct + index giraffe graph
# vg_giraffe_graph "" $ref chm13
# vg_giraffe_graph $variant_vcf $ref L1_vcfbub
# Construct + index vgmap graph
# bgzip Graph_genome_data/L1_annotation/L1_annotation.vcfbub.vcf
# tabix -p vcf Graph_genome_data/L1_annotation/L1_annotation.vcfbub.vcf.gz

# vg_map_graph "" $ref chm13
vg_map_graph $variant_vcf $ref L1_vcfbub

#K27ac
# graph_type=hprc-v1.1-mc-chm13
# pipeline=vg_giraffe
# fragment_length=120
# read_length=102
# unique_reads=13295105
# callpeaks "treatment_alignments.filtered.json" "" "${results_dir#$wd/}" "$graph_type" $pipeline $fragment_length $read_length $unique_reads


# # vg gamcompare -r 102 -T $wd/Graph_genome_chm13_trimmomatic/146_alignments.gam $wd/Graph_genome_vcfbub_trimmomatic/146_alignments.gam -t 64 > $wd/Graph_genome_data/146_alignments_compare.tsv
# # awk -F'\t' '$1 == 1 {correct_variant++} $1 == 0 {incorrect_variant++} END {
# #   print "Reads aligning better in the variant graph: " correct_variant;
# #   print "Reads aligning better in the linear graph: " incorrect_variant;
# # }' $wd/Graph_genome_data/146_alignments_compare.tsv


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
echo "Total runtime: ${elapsed} hours" >> slurm-${SLURM_JOB_ID}.err