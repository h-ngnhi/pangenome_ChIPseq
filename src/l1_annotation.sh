#!/bin/bash 
#SBATCH --job-name=vcfbub_L1_annotation
#SBATCH --output=slurm-%j.out  # %j will be replaced with the job ID
#SBATCH --error=slurm-%j.err   # %j will be replaced with the job ID
#SBATCH --account=def-bourqueg
#SBATCH --time=24:00:00
#SBATCH --mem=249G
#SBATCH --cpus-per-task=64

# Enable xtrace to log the commands to stderr
set -x

export wd=/home/hoangnhi/projects/def-bourqueg/hoangnhi/ZNF146-507-Analysis-on-Pangenome
export graph_type=Graph_genome_chm13_trimmomatic
graph_dir=$wd/$graph_type
graph_data_dir=$wd/Graph_genome_data
# cd $graph_dir

#1. Deconstruct the graph
# vg deconstruct -t 64 -P chr $graph_dir/hprc-v1.1-mc-chm13.gbz > $graph_dir/giraffe.vcf

#2. Edit the VCF file
# Define the input and output files
# awk 'BEGIN { FS=OFS="\t" }
# {
#     # Skip header lines
#     if ($0 ~ /^##/ || $0 ~ /^#CHROM/) {
#         print $0
#     } else {
#         # Replace ">" with "-" in the ID column (3rd column)
#         gsub(">", "-", $3)
#         print $0
#     }
# }' $graph_dir/giraffe.vcf > $graph_dir/giraffe_edited.vcf

#3. Filter top-level variant through vcfbub (No need to if already build the graph with vcfbub vcf file before)
# vcfbub --input $graph_dir/giraffe_edited.vcf --max-level 0 --max-ref-length 10000  > $graph_dir/giraffe_filt.vcf

#4. Make the graph biallelic because GraffiTE can't work with multiallelic variants
# module load bcftools/1.19
# zcat Graph_genome_full/hprc-v1.1-mc-chm13.raw.vcf.gz | bcftools norm -m- -o $graph_dir/giraffe_norm.vcf

#5. Filter variants that's less than 60bp
# awk 'BEGIN {FS="\t"; OFS="\t"} /^#/ {print $0; next} {if (length($4) >= 60 || length($5) >= 60) print $0}' $graph_dir/giraffe_norm.vcf > $graph_dir/giraffe_filt.vcf

#6. Check if variants are unique, rename if not
# awk '!/^#/ {print $3}' $graph_dir/giraffe_filt.vcf | sort | uniq -d
# Use awk to process the VCF file and rename the variant IDs
# awk 'BEGIN { OFS="\t"; id=1; } {if ($0 ~ /^#/) {print $0;} else {$3 = "var" id++;print $0;}}' $graph_dir/giraffe_filt.vcf > $graph_dir/giraffe_renamed.vcf


#5. Run GraffiTE

# module load nextflow/23.10.0
# module load apptainer
# nextflow run $wd/GraffiTE/main.nf \
#    --vcf $graph_dir/giraffe_renamed.vcf \
#    --TE_library $wd/GraffiTE_testset/human_DFAM3.6.fasta \
#    --reference $graph_data_dir/chm13.draft_v1.1.fasta \
#    --graph_method giraffe \
#    --mammal --cores 64 \
#    --repeatmasker_memory 249G \
#    --repeatmasker_time 48h \
#    -with-singularity $wd/graffite_latest.sif \
#    -profile cluster \
#    --genotype false

# #6. Filter L1 variants
# graph_dir=$wd/Graph_genome_vcfbub/
# bcftools view -i 'INFO/matching_classes ~ "LINE/L1"' $graph_dir/out/2_Repeat_Filtering/genotypes_repmasked_filtered.vcf -o $graph_dir/L1_annotation.vcf
# bcftools view -H $graph_dir/L1_annotation.vcf | wc -l
# # bcftools view -H $graph_dir/out/2_Repeat_Filtering/genotypes_repmasked_filtered.vcf | wc -l

# #7. Convert the graph to gfa file (don't work because odgi inject took really long time)
# vg convert -f $graph_dir/giraffe_index.gbz > $graph_dir/chm13v1graph.gfa
# odgi build -g $graph_dir/chm13v1graph.gfa -o $graph_dir/chm13v1graph.og
# linear_L1=$wd/Graph_genome_data/t2t.L1.noY.bed
# odgi inject -t 64 -i $graph_dir/chm13v1graph.og -b $linear_L1 -o $graph_dir/chm13v1graph_with_L1.og
# odgi paths -t 64 -i $graph_dir/chm13v1graph_with_L1.og --haplotypes -L > $graph_dir/touched_notes.txt

#8. Try the newest vg annotation (vg 1.61.0) (too long also :(((( )
chm13_L1=$wd/Graph_genome_data/t2t.L1.noY.bed
# vg annotate -b $chm13_L1 -x Graph_genome_vcfbub_trimmomatic/giraffe_index.gbz -F > $wd/Graph_genome_data/RMbed_L1_annotation.gaf
# vg deconstruct -t 64 -P chr $wd/Graph_genome_data/RMbed_L1_annotation.vg > $wd/Graph_genome_data/RMbed_L1_annotation.vcf

#9. Try vg find to extract node IDs associated with each LINE1 region
while read -r chrom start end l1_subfamily; do
    # Extract the subgraph for each region and convert to JSON format
    vg find -x Graph_genome_vcfbub_trimmomatic/giraffe_index.gbz -p "${chrom}:${start}-${end}" | \
    vg view -j - | \
    jq -r --arg l1_subfamily "$l1_subfamily" '
        .path[0].mapping as $mappings |
        [ 
            .node[] | 
            { id: .id, sequence: .sequence } 
        ] as $nodes |
        (
            # Generate the node traversal column
            ($mappings | 
            map(
                (if .position.is_reverse == true then "<" else ">" end) + (.position.node_id | tostring)
            ) | join(">")),

            # Generate the sequence column
            ($mappings | 
            map(
                ($nodes[] | select(.id == .position.node_id).sequence)
            ) | join(">")),

            # Include the L1 subfamily
            $l1_subfamily
        ) | @tsv
    '
done < $chm13_L1 > $wd/Graph_genome_data/RMbed_L1_vcf_node_traversals.txt



#10 Try vg paths 
# while read -r chrom start end family; do
#     echo "${chrom}:${start}-${end}#${family}"
# done < $chm13_L1 > $wd/Graph_genome_data/RMbed_L1_forvgpaths.txt

# while read -r region; do
#     chrom_region=$(echo "$region" | cut -d '#' -f 1)
#     family=$(echo "$region" | cut -d '#' -f 2)
    
#     # Append the family name as a GAF comment
#     vg paths -x Graph_genome_vcfbub_trimmomatic/giraffe_index.gbz -A -p "$chrom_region" | sed "s/$/\t#family:${family}/" >> $wd/Graph_genome_data/RMbed_L1_vcfbub_node_ids.gaf
# done < $wd/Graph_genome_data/RMbed_L1_forvgpaths.txt

########## The code below doens't work because i should not just concatenate the L1 variants from the CHM13 and the L1 variants from the graph genome

# #7. Convert the VCF file to BED file
# bcftools query -f '%CHROM\t%POS0\t%END\n' $graph_dir/L1_annotation.vcf > $graph_dir/L1_variants.bed
# bcftools query -f '%CHROM\t%POS\t%ALT\t%INFO/repeat_ids\n' L1_annotation.vcf | \
# awk -F'\t' 'BEGIN {OFS="\t"} {
#     split($4, ids, ",");
#     l1_ids = "";
#     for (i in ids) {
#         if (ids[i] ~ /^L1/) {
#             if (l1_ids != "") {
#                 l1_ids = l1_ids "," ids[i];
#             } else {
#                 l1_ids = ids[i];
#             }
#         }
#     }
#     if (l1_ids != "") {
#         print $1, $2-1, $2+length($3)-1, l1_ids;
#     }
# }' > L1_variants_type.bed

# #8. Combine L1 polymorphic and CHM13 L1 annotation
# chm13_L1=$wd/Graph_genome_data/t2t.L1.noY.bed
# cat $chm13_L1 L1_variants_type.bed | sort -k1,1 -k2,2n > L1_combined.bed

cat t2t.L1.noY.bed L1_variants.bed | sort -k1,1 -k2,2n > L1_combined.bed

################
# Read stats
# Remove chrY from the bed file
# awk '$1 != "chrY"' Genome/t2t.L1.bed > $graph_data_dir/t2t.L1.noY.bed

# # Extract MAPQ scores
# export graph_type=Graph_genome_vcfbub_trimmomatic
# graph_dir=$wd/$graph_type
# BAM_FILE=$graph_dir/146_alignments.bam
# BED_FILE=$graph_data_dir/t2t.L1.noY.bed
# INTERSECTED_BAM="$graph_dir/146_L1_alignments.rm.bam"
# MAPQ_SCORES="$graph_dir/mapq_scores_L1.rm.txt"
# echo $BAM_FILE
# echo $BED_FILE
# echo $INTERSECTED_BAM
# echo $MAPQ_SCORES

# #Intersect the alignment bam file with the L1 variants
# module load bedtools/2.31.0
# bedtools intersect -abam $BAM_FILE -b $BED_FILE > $INTERSECTED_BAM
# #Get the MAPQ and plot histogram:
# module load samtools/1.20
# samtools view $INTERSECTED_BAM | awk '{print $5}' > $MAPQ_SCORES

# #for all
# MAPQ_SCORES="$graph_dir/mapq_scores.txt"
# samtools view $BAM_FILE | awk '{print $5}' > $MAPQ_SCORES


# #Table of MAPQ scores
# total_reads=$(wc -l < "$MAPQ_SCORES")
# mapq_1=$(awk '$1 <= 5 {count++} END {print count+0}' "$MAPQ_SCORES")
# mapq_between=$(awk '$1 > 5 && $1 < 60 {count++} END {print count+0}' "$MAPQ_SCORES")
# mapq_60=$(awk '$1 == 60 {count++} END {print count+0}' "$MAPQ_SCORES")
# # Calculate percentages
# perc_mapq_1=$(echo "scale=2; ($mapq_1/$total_reads)*100" | bc)
# perc_mapq_between=$(echo "scale=2; ($mapq_between/$total_reads)*100" | bc)
# perc_mapq_60=$(echo "scale=2; ($mapq_60/$total_reads)*100" | bc)
# # Print the results
# echo -e "Category\tNumber of Reads\tPercentage"
# echo -e "mapq=1\t\t$mapq_1\t\t\t$perc_mapq_1%"
# echo -e "1<mapq<60\t$mapq_between\t\t\t$perc_mapq_between%"
# echo -e "mapq=60\t\t$mapq_60\t\t\t$perc_mapq_60%"
# echo -e "Total\t\t$total_reads\t\t\t100%"



# ##############
# # Peaks stats
# peak_file=/home/hoangnhi/projects/def-bourqueg/hoangnhi/ZNF146-507-Analysis-on-Pangenome/Graph_genome_chm13/callpeaks/graph_chm13_peaks.narrowPeak
# # awk '$1 != "chrY"' $peak_file > $graph_data_dir/146GenPipes.peaks.noY.bed
# bedtools intersect -a $peak_file -b $graph_data_dir/t2t.L1.noY.bed -u > Graph_genome_chm13/callpeaks/146peaks_L1.bed

# ref_l1="$graph_data_dir/t2t.L1.noY.bed" # repeatmasker file
# rm=$(wc -l $ref_l1) # total number of L1 peaks in the repeatmasker file
# peaks=$(wc -l "$graph_data_dir/146GenPipes.peaks.noY.bed")
# blacklist=Genome/Blacklist/t2t.excluderanges.bed

# # count overlap with L1 and calculate p_obs
# obs=$(wc -l "$graph_data_dir/146peaks_L1.bed") # total number of peaks overlapping with L1
# p_obs=$(echo ${obs%% *} / ${rm%% *} | bc -l) # proportion of peaks overlapping with L1
# ref=t2t
# p_shuffle=0 
# # shuffle 10 times, exclude blacklist region (not so significant but the paper method said so)
# for _ in {1..10}
# do 
#     bedtools shuffle -i "$graph_data_dir/146GenPipes.peaks.noY.bed" -g "Genome/human.$ref._noCHR.genome" -noOverlapping -maxTries 1000 -excl $blacklist > "$graph_data_dir/146GenPipes.peaks.noY.shuffle.bed"
#     awk '{print "chr"$0}' "$graph_data_dir/146GenPipes.peaks.noY.shuffle.bed" > temp_file && mv temp_file "$graph_data_dir/146GenPipes.peaks.noY.shuffle.bed"
#     # awk -i inplace '{print "chr"$0}' "$file.shuffle.bed" # output from bedtools doesn't have "chr"$0
#     sort -k1,1 -k2,2n "$graph_data_dir/146GenPipes.peaks.noY.shuffle.bed" -o "$graph_data_dir/146GenPipes.peaks.noY.shuffle.bed"
#     bedtools intersect -a "$graph_data_dir/146GenPipes.peaks.noY.shuffle.bed" -b $ref_l1 -u > "$graph_data_dir/146GenPipes.peaks.noY.shuffle.L1.bed"
#     shuffle=$(wc -l "$graph_data_dir/146GenPipes.peaks.noY.shuffle.L1.bed") # total number of peaks overlapping with L1 after shuffling
#     echo "${shuffle%% *}"
#     p_shuffle=$(echo $p_shuffle + ${shuffle%% *} / ${rm%% *} | bc -l) # summary of proportion of peaks overlapping with L1 after shuffling
#     echo "$p_shuffle"
# done
# p_shuffle=$(echo $p_shuffle/10 | bc -l) # average of proportion of peaks overlapping with L1 after shuffling
# enr=$(echo $p_obs / $p_shuffle | bc -l)
# echo "enrichment $enr"