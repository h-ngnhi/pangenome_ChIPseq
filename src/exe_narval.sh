#!/bin/bash 

export wd=$(pwd)
# Define parameters
export mark_triple=("K27_FLU 17 7") # e.g., K27_FLU, 146, 507, or iPSC_K27 then k and w if test for giraffe param Eg. "146 17 7"
pipeline=("vg_giraffe")       # vg_giraffe or vg_map
ref=("vcfbub")               # e.g., chm13, L1_vcfbub, vcfbub, or hprc-v1.1-mc-chm13
steps="4 5 6"             # Steps of the pipeline
        # Steps of the pipeline
            # 1. Construct the graph
            # 2. Split the graph
            # 3. Align the reads
            # 4. Calculate parameters for peak calling
            # 5. Split json to chromosomes
            # 6. Call peaks + Callpeaks_whole_genome_from_p_values + Find linear
            # 7. Callpeaks_whole_genome_from_p_values + Find linear
      
mapq_troubleshoot=("" "mapq60")        # mapq60 or inject or hard_hit_cap or "" - for troubleshooting mapq
        # mapq60 is to change all alignments with MAPQ < 60 to 60 then call peaks
        # inject is to inject the linear reads into the graph then call peaks
        # hard_hit_cap is to increase the hard hit reads to be larger than default 500

mapq_filter=("10")      #"30" at default will remove all alignments with MAPQ < 30

# # Export variables for the job (they will be available in the sbatch command)
export pipeline ref steps mapq_troubleshoot

# Name of the log file
LOG_FILE="$wd/pangenome_ChIPseq/april_june_2025.txt"
TEMP_LOG="$wd/pangenome_ChIPseq/temp_job_log.txt"

echo "=== $(date '+%Y-%m-%d %H:%M:%S') ===" > "$TEMP_LOG"
# Use GNU Parallel to submit up to 5 jobs concurrently.
# Use double quotes in the parallel command so variables get expanded.
parallel -j5 --linebuffer "
    IFS=' ' read mark k w <<< {1};
    JOB_OUT=\$(sbatch -J chipseq_\${mark}_{2}_{3} \$wd/pangenome_ChIPseq/src/exe.sh \$mark {2} {3} \"\$steps\" {4} {5} \$k \$w)
    JOB_ID=\$(echo \"\$JOB_OUT\" | awk '{print \$4}')
    echo \"\$mark {2} {3} \$steps {4} mapq{5} k=\$k w=\$w -> job_id:\$JOB_ID\"
" ::: "${mark_triple[@]}" ::: "${pipeline[@]}" ::: "${ref[@]}" ::: "${mapq_troubleshoot[@]}" ::: "${mapq_filter[@]}" >> "$TEMP_LOG"

# parallel -j5 --linebuffer "
#     JOB_OUT=\$(sbatch -J bwa_{1} \$wd/temp2.sh {1})
#     JOB_ID=\$(echo \"\$JOB_OUT\" | awk '{print \$4}')
#     echo \"\$(date '+%Y-%m-%d %H:%M:%S') {1} -> job_id:\$JOB_ID\"
# " ::: "${markname[@]}" >> "$TEMP_LOG"

# Now, prepend the TEMP_LOG content at the beginning of LOG_FILE.
tmp=$(mktemp)
cat "$TEMP_LOG" > "$tmp"
cat "$LOG_FILE" >> "$tmp"
mv "$tmp" "$LOG_FILE"
rm "$TEMP_LOG"


# vg find -x Pangenomes/vg_giraffe/vcfbub/vcfbub.gbz -p chr13:1000-2500 -E > Pangenomes/vg_giraffe/vcfbub/giraffe_vcfbub_chr13_1500_2500.vg
# vg convert -f Pangenomes/vg_giraffe/vcfbub/giraffe_vcfbub_chr13_1500_2500.vg > Pangenomes/vg_giraffe/vcfbub/giraffe_vcfbub_chr13_1500_2500.gfa
# vg find -x Pangenomes/vg_giraffe/vcfbub/vcfbub.gbz -p chr12:1500-2500 -E

# Access MAPQ GenPipes
# chipseq_gp K27_FLU $wd/results/K27_FLU "-f"
# prim_bam results/K27_FLU/linear_chm13/alignment/H3K27AC_treatment1/H3K27AC/H3K27AC_treatment1.H3K27AC.sorted.dup.bam
# align_stats results/K27_FLU/linear_chm13/alignment/H3K27AC_treatment1/H3K27AC/H3K27AC_treatment1.H3K27AC.sorted.dup.filtered.cleaned.bam



# enrichment_calc() {
#     local bed_path=$1   # to locate the bed file
#     local bed_file=$2   # to take the file name and create downstream files
#     local csv_file=$3   # to write the enrichment csv file
#     local te=$4        # to use which ref genome
#     local blacklist=$5  # to use which blacklist file
#     local row_title=$6  # the title of each calculation (1st column)
    
#     module load bedtools
#     ref=t2t
#     ref_l1="$wd/Genome/t2t.$te.bed" # repeatmasker file
#     rm=$(wc -l $ref_l1) # total number of L1 peaks in the repeatmasker file
#     file=$bed_path/$bed_file
#     echo $file
#     peaks=$(wc -l "$file.narrowPeak.bed")
#     sort -k1,1 -k2,2n "$file.narrowPeak.bed" > "$file.sorted.bed"
#     if ! grep -Eq '^(chr)' "$file.sorted.bed"; then
#         awk '{print "chr"$0}' "$file.sorted.bed" > temp_file && mv temp_file "$file.sorted.bed"
#         # awk -i inplace '{print "chr"$0}' "$file.sorted.bed" # apply only to hg19 GenPipes file because somehow they don't includr "chr"
#     fi

#     # count overlap with L1 and calculate p_obs
#     # NEED TO OPTIMIZE: I should've count the peaks only instead of creating a new bed file
#     bedtools intersect -a $ref_l1 -b "$file.sorted.bed" -u > "$bed_path/all_intersection.$ref.$te.bed"
#     obs=$(wc -l "$bed_path/all_intersection.$ref.bed") # total number of peaks overlapping with L1
#     p_obs=$(echo ${obs%% *} / ${rm%% *} | bc -l) # proportion of peaks overlapping with L1

#     p_shuffle=0 
#     # shuffle 10 times, exclude blacklist region (not so significant but the paper method said so)
#     for _ in {1..10}
#     do 
#         bedtools shuffle -i "$file.sorted.bed" -g "Genome/human.$ref._noCHR.genome" -noOverlapping -maxTries 1000 -excl $blacklist > "$file.shuffle.bed"
#         awk '{print "chr"$0}' "$file.shuffle.bed" > temp_file && mv temp_file "$file.shuffle.bed"
#         # awk -i inplace '{print "chr"$0}' "$file.shuffle.bed" # output from bedtools doesn't have "chr"$0
#         sort -k1,1 -k2,2n "$file.shuffle.bed" -o "$file.shuffle.bed"
#         bedtools intersect -a "$file.shuffle.bed" -b $ref_l1 -u > "$bed_path/all_intersection_shuffle.$ref.bed"
#         shuffle=$(wc -l "$bed_path/all_intersection_shuffle.$ref.bed") # total number of peaks overlapping with L1 after shuffling
#         echo "${shuffle%% *}"
#         p_shuffle=$(echo $p_shuffle + ${shuffle%% *} / ${rm%% *} | bc -l) # summary of proportion of peaks overlapping with L1 after shuffling
#         echo "$p_shuffle"
#     done

#     p_shuffle=$(echo $p_shuffle/10 | bc -l) # average of proportion of peaks overlapping with L1 after shuffling

#     enr=$(echo $p_obs / $p_shuffle | bc -l)
#     echo "enrichment $enr"
    
#     echo -e "${row_title}_${gene}_${ref}\t${peaks%% *}\t${obs%% *}\t$enr" >> "$csv_file"
# }
# bed_path=results/iPSC_K27/linear_chm13/peak_call/ipsc_K27/ipsc_K27/   # to locate the bed file
# bed_file=ipsc_K27.ipsc_K27_peaks   # to take the file name and create downstream files
# te=sva       # to use which ref genome
# csv_file=results/iPSC_K27/linear_chm13/${te}_enr.csv   # to write the enrichment csv file
# blacklist=Genome/Blacklist/t2t.excluderanges.bed  # to use which blacklist file
# row_title=iPSC_${te}_linear_chm13   # to use which row title
# enrichment_calc $bed_path $bed_file $csv_file $te $blacklist $row_title
# bed=results/iPSC_K27/linear_chm13/peak_call/ipsc_K27/ipsc_K27/all_intersection.t2t.$te.bed
# wc -l $bed
# awk '{print $4}' $bed | sort | uniq -c | sort -nr | head -n 10

# java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 4 \
#   results/K27_FLU/trim/H3K27AC_treatment1/H3K27AC/Treatment1.trim.pair1.fastq.gz results/K27_FLU/trim/H3K27AC_treatment1/H3K27AC/Treatment1.trim.pair2.fastq.gz \
#   results/K27_FLU/trim/H3K27AC_treatment1/H3K27AC/Treatment1_15.trim.pair1.fastq.gz results/K27_FLU/trim/H3K27AC_treatment1/H3K27AC/Treatment1_15.trim.unpair1.fastq.gz \
#   results/K27_FLU/trim/H3K27AC_treatment1/H3K27AC/Treatment1_15.trim.pair2.fastq.gz results/K27_FLU/trim/H3K27AC_treatment1/H3K27AC/Treatment1_15.trim.unpair2.fastq.gz \
#   HEADCROP:15
# java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 4 \
#   data/iPSC/p_N_iPSC_K27ac_r2_S20_R1_001.fastq.gz data/iPSC/p_N_iPSC_K27ac_r2_S20_R2_001.fastq.gz \
#   data/iPSC/Treatment2.trim.pair1.fastq.gz data/iPSC/Treatment2.trim.unpair1.fastq.gz \
#   data/iPSC/Treatment2.trim.pair2.fastq.gz data/iPSC/Treatment2.trim.unpair2.fastq.gz \
#   HEADCROP:8
# java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 4 \
#   data/iPSC/p_N_iPSC_Input_r1_S18_R1_001.fastq.gz data/iPSC/p_N_iPSC_Input_r1_S18_R2_001.fastq.gz \
#   data/iPSC/Input.trim.pair1.fastq.gz data/iPSC/Input.trim.unpair1.fastq.gz \
#   data/iPSC/Input.trim.pair2.fastq.gz data/iPSC/Input.trim.unpair2.fastq.gz \
#   HEADCROP:8