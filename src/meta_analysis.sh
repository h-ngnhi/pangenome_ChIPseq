#!/bin/bash
export wd=/home/hoangnhi/projects/def-bourqueg/hoangnhi/ZNF146-507-Analysis-on-Pangenome
cd $wd
module load StdEnv/2020
module load samtools/1.10
module load bedtools/2.30.0
module load kentutils/401 # bigWigToBedGraph, bedGraphToBigWig
module load gcc/9.3.0 #required for meme
# module load openmpi/4.0.3   #required for meme
module load meme/5.5.5



##################################################
# Plot ENCODE data
intersect(){
    local file=${1%.bigWig}
    local ref=$2
    local ref_l1="$wd/Genome/$ref.L1.bed" 

    bigWigToBedGraph "$file.bigWig" "$file.bedGraph"
    # sort -k1,1 -k2,2n "$file.bedGraph" > "$file.sorted.bedGraph" #sorting always get killed
    bedtools intersect -a "$file.bedGraph" -b $ref_l1 -u > "$file.L1.bedGraph"
    bedGraphToBigWig "$file.L1.bedGraph" "$wd/meta_analysis_files/$ref.chrom.sizes" "$file.L1.bigWig"
}

heatmap_creator(){
    local bw=$1 #BigWig file
    local name=$2 #Name of the output file
    local disrupted=$3 #Disrupted LINE-1 file
    local intact=$4 #Intact LINE-1 file
    local color=$5 #Color of the heatmap
    local title=$6 #Title of the heatmap
    local height=$7 #Height of the heatmap
    local width=$8 #Width of the heatmap
    # local color=$5 #Color of the heatmap
    
    computeMatrix reference-point --referencePoint center \
              -b 4000 -a 4000 \
              -R $disrupted $intact \
              -S $bw \
              --binSize 10 \
              --missingDataAsZero \
              -o "$meta_directory/matrix.gz"
    plotHeatmap -m "$meta_directory/matrix.gz" \
            -out $name \
            --colorMap $color \
            --regionsLabel "Disrupted LINE-1" "Intact LINE-1" \
            --sortRegions descend \
            --plotTitle "$title" \
            --heatmapHeight $height --heatmapWidth $width \
            --refPointLabel "LINE-1 Center"
}


# tail -n +2 meta_analysis_files/disrupted.hg38.bed > meta_analysis_files/disrupted_no_header.hg38.bed
# tail -n +2 meta_analysis_files/intact.hg38.bed > meta_analysis_files/intact_no_header.hg38.bed
# NOTE that disrupted and intact files are converted from hg18 to hg19 using UCSC liftOver website tool
date=$(date "+%Y-%m-%d")


###
for ref in hg19 hg38; do
    disrupted=meta_analysis_files/disrupted_no_header.$ref.bed
    intact=meta_analysis_files/intact_no_header.$ref.bed
    # Generate heatmap for 146 and 507
    for gene in 146 507; do
        meta_directory=$wd/meta_analysis_files/ENCODE/$gene
        mkdir -p $meta_directory

        # 146_unique
        # Intersect bigWig file with LINE-1 elements
        # bw_file=$(find $meta_directory -type f -name "*.${ref}.bigWig")
        # intersect $bw_file $ref
        # # Generate heatmap
        # heatmap_creator "${bw_file%.bigWig}.L1.bigWig" "$meta_directory/$date.unique_heatmap.$ref.png" $disrupted $intact Blues "Unique LINE-1 elements"

        # 146_primary (+multimap)
        bam_files=($(find $meta_directory -type f -name "*.${ref}.bam"))
        # Merge and sort bam files
        if [ ! -d "$meta_directory/$merged.$ref.sorted.bam" ]; then
            samtools merge -f $meta_directory/merged.$ref.bam ${bam_files[0]} ${bam_files[1]}
            samtools sort -o $meta_directory/merged.$ref.sorted.bam $meta_directory/merged.$ref.bam
            # samtools index $meta_directory/merged.$ref.sorted.bam # Index is not necessary for bedtools
        fi
        echo "Merged and sorted"
        # Intersect bam file with LINE-1 elements
        bedtools intersect -a "$meta_directory/merged.$ref.sorted.bam" -b "$wd/Genome/$ref.L1.bed"  -u > "$meta_directory/merged.$ref.sorted.L1.bam"
        # Index the intercepted bam file
        samtools index $meta_directory/merged.$ref.sorted.L1.bam
        # Generate coverage file
        # "coverage files of fold change over control coverage were also generated without filtering for uniqueness 
        # using deepTools bamCoverage and the –normalizeTo1X parameters"
        bamCoverage --bam $meta_directory/merged.$ref.sorted.L1.bam --outFileName $meta_directory/merged_bam_normalized.$ref.bw \
                    --binSize 15 --effectiveGenomeSize 2700000000 --normalizeUsing RPGC --ignoreDuplicates --extendReads
        # Generate heatmap
        heatmap_creator "$meta_directory/merged_bam_normalized.$ref.bw" "$meta_directory/$date.primary_heatmap.$ref.png" $disrupted $intact Blues "Primary (+multimap) LINE-1 elements" 5 20
    done

    # Generate mappability heatmap
    map_file=$(find $wd/meta_analysis_files/ -maxdepth 1 -type f -name "*.${ref}.bigWig" )
    if [ ! -d "$meta_directory/${map_file%.bigWig}.L1.bigWig" ]; then
            intersect $map_file $ref
    fi
    heatmap_creator "${map_file%.bigWig}.L1.bigWig" "$wd/meta_analysis_files/$date.mappability_heatmap.$ref.png" $disrupted $intact Greys_r "Mappability" 5 20
done

##################################################
# Plot GENPIPES data - primary
type="design2"
ref=hg19
# gene=507
gene=146
map=unique
disrupted=meta_analysis_files/disrupted_no_header.$ref.bed
intact=meta_analysis_files/intact_no_header.$ref.bed


meta_directory=$wd/${gene}_${type}_${ref}/alignment/$gene/ZNF$gene
meta_plot_directory=$wd/meta_analysis_files/GenPipes/$gene
mkdir -p $meta_plot_directory

for ref in hg19 hg38 t2t; do
    disrupted=meta_analysis_files/disrupted_no_header.$ref.bed
    intact=meta_analysis_files/intact_no_header.$ref.bed
    # Generate heatmap for 146 and 507
    for gene in 146 507; do
        meta_directory=$wd/${gene}_${type}_${ref}/alignment/$gene/ZNF$gene
        meta_plot_directory=$wd/meta_analysis_files/GenPipes/$gene

        for map in unique primary; do
            module load StdEnv/2020
            module load samtools/1.10
            module load bedtools/2.30.0
            bam_file=$(get_file "GenPipes_bam_$map" $gene $ref)
            echo $bam_file
            # Need to add "chr" prefix to the chromosome names in the GenPipes bam files
            if [ $ref == hg19 ]; then
                # Extract the header of the BAM file
                echo $ref
                # Function to clean and format the chromosome names
                clean_header() {
                    awk '{
                        if ($1 ~ /^@SQ/) {
                            sub(/SN:(chr)+/, "SN:chr", $0);  # Remove any extra "chr" prefixes
                            sub(/SN:([0-9XYMT])/, "SN:chr\\1", $0);  # Add "chr" prefix if missing
                        }
                        print $0
                    }' "$1" > "$1.cleaned"
                    mv "$1.cleaned" "$1"
                }
                # Check if the header contains incorrectly formatted chromosome names and fix them
                clean_header "$meta_directory/header.sam"
                # Reheader the BAM file with the cleaned header
                samtools reheader "$meta_directory/header.sam" $bam_file > "$meta_directory/$gene.ZNF$gene.sorted.reheader.bam"
                # Replace the original BAM file with the reheadered BAM file
                mv "$meta_directory/$gene.ZNF$gene.sorted.reheader.bam" $bam_file                
                # Index the BAM file
                samtools index $bam_file
            fi
            L1_bed=$(get_file "L1_bed" $gene $ref)
            merged_bam=$meta_directory/merged.$ref.$map.sorted.L1.bam
            bedtools intersect -abam $bam_file -b $L1_bed -u > $merged_bam
            head $merged_bam
            # Index the intercepted bam file
            samtools index $merged_bam
            # Generate coverage file

            module purge 
            module load mugqic/deepTools/3.5.1 
            bw_file=$meta_directory/$ref.$map.L1.bw
            bamCoverage --bam $merged_bam --outFileName $bw_file \
                        --binSize 15 --effectiveGenomeSize 2700000000 --normalizeUsing RPGC --ignoreDuplicates --extendReads
            # version `XCRYPT_2.0' not found error 
            
            # Generate heatmap
            heatmap_creator "$meta_directory/merged_bam_normalized.$ref.bw" "$meta_plot_directory/$date.primary_heatmap.$ref.png" $disrupted $intact Blues "Primary (+multimap) LINE-1 elements"
        fi
    done
done


a=/home/hoangnhi/projects/def-bourqueg/hoangnhi/ZNF146-507-Analysis-on-Pangenome/146_design2_hg19/peak_call/146/ZNF146/146.ZNF146_peaks.sorted.bed
b=/home/hoangnhi/projects/def-bourqueg/hoangnhi/ZNF146-507-Analysis-on-Pangenome/meta_analysis_files/L1HS.hg19.bed
bedtools intersect -a $a -b $b  -u > "$meta_directory/hg19.GenPipes.L1HS.bed"

##################################################
### Analyse only unique T2T regions

a=/home/hoangnhi/projects/def-bourqueg/hoangnhi/ZNF146-507-Analysis-on-Pangenome/146_design2_hg19/peak_call/146/ZNF146/146.ZNF146_peaks.sorted.bed
b=/home/hoangnhi/projects/def-bourqueg/hoangnhi/ZNF146-507-Analysis-on-Pangenome/meta_analysis_files/L1HS.hg19.bed
bedtools intersect -a $a -b $b  -u > "$meta_directory/hg19.GenPipes.L1HS.bed"
            
##################################################
###Coverage of peaks aligned to L1HS (L1PA1)
type="design2"
ref=hg19
gene=146

# ENCODE
bed_file=$wd/GenPipes_sets/data/ENCODE/$gene/all_intersection.$ref.bed
meta_plot_directory=$wd/meta_analysis_files/ENCODE/$gene



# GenPipes
fasta_file=$wd/Genome/$ref.fa
bed_file=$wd/${gene}_${type}_${ref}/peak_call/$gene/ZNF$gene/all_intersection.$ref.bed
meta_plot_directory=$wd/meta_analysis_files/GenPipes/$gene


# Get fasta file
bedtools getfasta -fi $fasta_file -bed $bed_file -fo $meta_plot_directory/$gene.ZNF${gene}_peaks.fa

# Run meme to discover motifs
meme $meta_plot_directory/$gene.ZNF${gene}_peaks.fa -dna -mod zoops -nmotifs 1 -minw 7 -maxw 13 -maxsize 15000000 -oc $meta_plot_directory/meme_output
# -dna: Specifies that the input sequences are DNA.
# -mod zoops: Zero or One Occurrence Per Sequence, which means each motif is allowed to appear no more than once in any given sequence.
# -nmotifs 1: The number of different motifs to find.
# -minw 6 and -maxw 13: Minimum and maximum motif width to consider.
# -maxsize 1000000: The maximum total number of bases to use from the input file.

    # Run fimo to scan motifs
    fimo_output=$meta_plot_directory/$ref.fimo_output_L1
    mkdir -p $fimo_output
    fimo --text --thresh 1e-4 --oc $fimo_output $meme_output/meme.txt $peak_L1HS_fa > $fimo_output/fimo.txt
    # --text: Output the results in text format.
    # --thresh 1e-4: The threshold for reporting motif matches. This is the p-value threshold for reporting motif matches.
    # --oc fimo_output: The output directory for FIMO results.
    # $gene.ZNF${gene}_peaks.fa: The input sequences to scan for motifs.
    # meme_out/meme.txt: The motif file from MEME.



# less $wd/meta_analysis_files/ENCODE/507/hg19.ZNF507_peaks.L1HS.bed
genes=("146" "507")
refs=("t2t") # will modify  later to include "T2T"
type="design2" 
pipelines=("ENCODE" "GenPipes")
for gene in "${genes[@]}"; do
    for ref in "${refs[@]}"; do
        for pipeline in "${pipelines[@]}"; do
            echo ""
            echo "Gene: $gene, Reference: $ref, Pipeline: $pipeline"
            # ENCODE
            # motif_ana $gene $ref "ENCODE" 

            # GenPipes
            motif_ana $gene $ref "GenPipes"

            fimo_output=$wd/meta_analysis_files/$pipeline/$gene/$ref.fimo_output_L1
            motif_id=$(awk 'NR==2 {print $1}' "$fimo_output/fimo.txt")
            echo "Motif ID: $motif_id"
            # Read the consensus sequence from the second line of L1HS.fa
            consensus_sequence=$(sed -n '2p' "$(get_file "L1HS_fa")")
            # echo "$consensus_sequence" | awk -v motif="${motif_id,,}" '{ print index($0, motif) }'
            echo "$consensus_sequence" | awk -v motif="${motif_id,,}" \
              '{ idx = index($0, motif); if(idx > 0) { print idx, idx + 12 } else { print "Motif not found" } }'

        done
    done
done

#######################################################
### Calculate coverage of peaks aligned to L1HS (L1PA1) for plot
type="design2"
ref=hg19
gene=146
pipeline="GenPipes"

genes=("146" "507")
refs=("hg19" "hg38") # will modify  later to include "T2T"
type="design2" 
pipelines=("GenPipes")
for gene in "${genes[@]}"; do
    for ref in "${refs[@]}"; do
        for pipeline in "${pipelines[@]}"; do
            echo ""
            echo "Gene: $gene, Reference: $ref, Pipeline: $pipeline"

            meta_plot_directory=$wd/meta_analysis_files/$pipeline/$gene
            ref_fa=$wd/Genome/$ref.fa
            peak_L1HS_bed=$meta_plot_directory/$ref.ZNF${gene}_peaks.L1.bed
            peak_L1HS_fa=$meta_plot_directory/$ref.ZNF${gene}_peaks.L1.fa
            bedtools getfasta -fi $ref_fa -bed $peak_L1HS_bed -fo $peak_L1HS_fa -s

            # Align to Consensus sequence
            L1HS_fa=$(get_file "L1HS_fa")
            aligned_sam=$meta_plot_directory/$ref.aligned_reads.sam
            bwa mem -B 1 -O 1 -E 1 -T 1 -t 12 -k 10 $L1HS_fa $peak_L1HS_fa > $aligned_sam

            # Convert SAM to BAM and sort
            samtools view -Sb $aligned_sam | samtools sort -o $meta_plot_directory/$ref.aligned_reads.sorted.bam

            # Calculate coverage
            samtools depth -a $meta_plot_directory/$ref.aligned_reads.sorted.bam > $meta_plot_directory/$ref.L1_coverage.txt
            echo "Coverage file: $meta_plot_directory/$ref.L1_coverage.txt"
        done  
    done
done





# ##################################################
# # Use deepTools bamCoverage to generate coverage files
# input_bam=GenPipes_sets/data/146_eGFP_HEK293/ENCFF830EGQ_rep1.hg19.bed
# bamCoverage -b $input_bam -o $meta_directory/coverage_rep1.bw --normalizeUsing CPM 

# GenPipes_sets/data/146_eGFP_HEK293/ENCFF830EGQ_rep1.hg19.bed




# export wd=/home/hoangnhi/projects/def-bourqueg/hoangnhi/ZNF146-507-Analysis-on-Pangenome/146_design_hg19
# # To load mugqic module
# # export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
# # module use $MUGQIC_INSTALL_HOME/modulefiles
# # List of modules are in chipseq.base.ini https://bitbucket.org/mugqic/genpipes/src/ffdc708fdc0e8d3585dff0acbbc1dea2b9b5b723/pipelines/chipseq/chipseq.base.ini?at=master#chipseq.base.ini-36:55,64,69,153:155
# # module load mugqic/MACS2/2.2.7.1
# # module load mugqic/python/3.10.4

# # Starting point is the result BAM files from SAMbamba View Filter step from GenPipes
# # 146_design_hg19/alignment/146Rep1/ZNF146/146Rep1.ZNF146.sorted.dup.filtered.cleaned.bam is output after bedtools blacklist filter
# cd $wd
# chip_bam=alignment/146Rep1/ZNF146/146Rep1.ZNF146.sorted.dup.filtered.cleaned.bam
# input_bam=alignment/146Input/Input/146Input.Input.sorted.dup.filtered.cleaned.bam
# macs2 callpeak -t $chip_bam -c $input_bam -f BAM -g hs -n 146Rep1.ZNF146_treat_pileup.bdg --outdir macs2_output -B -q 0.01
type="design2"
ref=hg19
gene=146

### HEATMAP CREATION
### Create fold change BigWig files for GenPipes ONLY (ENCODE already has fold change BigWig files)
# Code is inspired from ENCODE
# 1. Convert BAM files to BedGraph
meta_plot_directory=$wd/meta_analysis_files/GenPipes/$gene
pipeline_bam=$(get_file "GenPipes_bam_unique" $gene $ref)
pipeline_bdg=$meta_plot_directory/$gene.$ref.bdg
bamCoverage --bam $pipeline_bam --outFileName $pipeline_bdg --binSize 10 --outFileFormat bedgraph

ctrl_bam=$(get_file "GenPipes_bam_ctrl" $gene $ref)
ctrl_bdg=$meta_plot_directory/Input.$ref.bdg
bamCoverage --bam $ctrl_bam --outFileName $ctrl_bdg --binSize 10 --outFileFormat bedgraph

# 2. Calculate Fold Enrichment using MACS2
macs2 bdgcmp -t $pipeline_bdg -c $ctrl_bdg --outdir $meta_plot_directory -o $ref.unique_FE.bdg -m FE

# 3. Correcting Coordinates and Formatting
CHRSIZEFILE=$(get_file "chrsize_$ref") #<path_of_file_containing_chromosome_sizes>
slopBed -i $meta_plot_directory/$ref.unique_FE.bdg -g ${CHRSIZEFILE} -b 0 | awk '{if ($3 != -1) print $0}' \
  | bedClip stdin ${CHRSIZEFILE} $meta_plot_directory/$ref.unique.fc.signal.bedgraph

# 4. Cleanup intermediate files
rm -f $meta_plot_directory/$ref.unique_FE.bdg

# 5. Convert BedGraph to BigWig
bedGraphToBigWig $meta_plot_directory/$ref.unique.fc.signal.bedgraph ${CHRSIZEFILE} $meta_plot_directory/$ref.unique.fc.signal.bw

# bamCoverage -b ChIP_filtered.bam -o ChIP.bw --normalizeUsing RPKM --ignoreDuplicates

# # Control sample
# bamCoverage -b Input_filtered.bam -o Input.bw --normalizeUsing RPKM --ignoreDuplicates

############
#Plot histogram of mapping quality scores

# Define your input files
BAM_FILE=Graph_genome_full/146_alignments.bam
BED_FILE=Genome/t2t.L1.bed
INTERSECTED_BAM="Graph_genome/146_L1_alignments.bam"
MAPQ_SCORES="Graph_genome/mapq_scores.txt"

# Step 1: Intersect BAM with BED
bedtools intersect -a $BAM_FILE -b $BED_FILE -u > $INTERSECTED_BAM

# Step 2: Extract MAPQ Scores
samtools view $INTERSECTED_BAM | awk '{print $5}' > $MAPQ_SCORES

# Step 3: Plot the histogram using Python
python3 - <<EOF
import pandas as pd
import matplotlib.pyplot as plt

# Load the mapping quality scores
mapq_scores = pd.read_csv('$MAPQ_SCORES', header=None, names=['MAPQ'])

# Plot the histogram
plt.hist(mapq_scores['MAPQ'], bins=50, edgecolor='black')
plt.title('Distribution of Mapping Quality Scores')
plt.xlabel('MAPQ Score')
plt.ylabel('Frequency')
plt.grid(True)
plt.show()
EOF
