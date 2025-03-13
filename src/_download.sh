#!/bin/bash
#Scripts to download ZNF146 and ZNF507 raw sequencing data from ENCODE project
export wd=/home/hoangnhi/projects/def-bourqueg/hoangnhi/ZNF146-507-Analysis-on-Pangenome
cd $wd
download_directory=GenPipes_sets/data

download_read() {
  echo "#to be edit later"
}

download_published_hprc() {
  cd Pangenome/vg_giraffe/hprc-v1.1-mc-chm13 || exit
  wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-chm13/hprc-v1.1-mc-chm13.gbz
  wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-chm13/hprc-v1.1-mc-chm13.dist
  wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-chm13/hprc-v1.1-mc-chm13.min
}

download_ENCODE() {
  local experiment_id=$1
  local folder=$2

  file_info=$(curl -s "https://www.encodeproject.org/experiments/${experiment_id}/?format=json" | jq -c '
                  .files[] | 
                  select(.file_type == "bed narrowPeak" and .output_type == "optimal IDR thresholded peaks") | 
                  {file_url: .href, genome_assembly: .assembly}')
  while read -r line; do
    file_url=$(echo "$line" | jq -r '.file_url')
    genome_assembly=$(echo "$line" | jq -r '.genome_assembly')
    file_name="${download_directory}/${folder}/${genome_assembly}.$(basename "$file_url")"
    curl -o "$file_name" -J -L "https://www.encodeproject.org${file_url}"
    gunzip -c "$file_name" > "${file_name%.bed.gz}.narrowPeak.bed" && rm "$file_name"
  done <<< "$file_info" 
}

#ZNF146/ONF-eGFP-HEK293
#To be stored in data/146_eGFP_HEK293 folder
wget https://www.encodeproject.org/files/ENCFF441XRJ/@@download/ENCFF441XRJ.fastq.gz #Replicate 1; Read 1
wget https://www.encodeproject.org/files/ENCFF382JJC/@@download/ENCFF382JJC.fastq.gz #Replicate 1; Read 2
wget https://www.encodeproject.org/files/ENCFF209LAN/@@download/ENCFF209LAN.fastq.gz #Replicate 2; Read 1
wget https://www.encodeproject.org/files/ENCFF027IUI/@@download/ENCFF027IUI.fastq.gz #Replicate 2; Read 2

# Download 24 controls 
ENCODE_MATCHED_SET_ID="ENCSR598OXI"
# Fetch the list of experiments in the matched set           
experiment_ids=$(curl -s "https://www.encodeproject.org/matched-sets/$ENCODE_MATCHED_SET_ID/?format=json" | \
                 jq --arg ENCODE_MATCHED_SET_ID "$ENCODE_MATCHED_SET_ID" \
                 '.. | strings | select(startswith("ENCSR") and . != $ENCODE_MATCHED_SET_ID)' | \
                 grep -o 'ENCSR[^ "]*' | sort -u)

for experiment_id in $experiment_ids; do
  # Fetch the list of FASTQ file URLs for the experiment
  fastq_urls=$(curl -s "https://www.encodeproject.org/experiments/$experiment_id/?format=json" | jq -r '.files[] | select(.output_type == "reads") | .href')
  # Download the FASTQ files for the experiment
    i=0
    for fastq_url in $fastq_urls; do
        ((i++))
        wget -O "data/146_eGFP_HEK293/controls/Input$i.$experiment_id.${fastq_url: -20}" "https://www.encodeproject.org$fastq_url"
    done
done

# Download ENCODE processed data
folder=ENCODE/ZNF146
mkdir -p "$download_directory/$folder" #in case the folder is not created yet
download_ENCODE ENCSR689YFA $folder


#ZNF507-eGFP-K562 (Checked to be same as paper)
#To be stored in data/507_eGFP_K562 folder
wget https://www.encodeproject.org/files/ENCFF306MYT/@@download/ENCFF306MYT.fastq.gz #Replicate 1; Read 1
wget https://www.encodeproject.org/files/ENCFF574SDW/@@download/ENCFF574SDW.fastq.gz #Replicate 1; Read 2
wget https://www.encodeproject.org/files/ENCFF240EGC/@@download/ENCFF240EGC.fastq.gz #Replicate 2; Read 1
wget https://www.encodeproject.org/files/ENCFF905PWR/@@download/ENCFF905PWR.fastq.gz #Replicate 2; Read 2
#Download control ENCSR068GPE
wget https://www.encodeproject.org/files/ENCFF445KUN/@@download/ENCFF445KUN.fastq.gz #Read 1
wget https://www.encodeproject.org/files/ENCFF091IYB/@@download/ENCFF091IYB.fastq.gz #Read 2
#Download ENCODE processed data file
folder=ENCODE/ZNF507
mkdir -p "$download_directory/$folder" #in case the folder is not created yet
download_ENCODE ENCSR598TIR $folder

# experiment_id="ENCSR598TIR"
# file_info=$(curl -s "https://www.encodeproject.org/experiments/${experiment_id}/?format=json" | jq -c '
#                   .files[] | 
#                   select(.file_type == "bed narrowPeak" and .output_type == "optimal IDR thresholded peaks") | 
#                   {file_url: .href, genome_assembly: .assembly}')
# while read -r line; do
#   file_url=$(echo "$line" | jq -r '.file_url')
#   genome_assembly=$(echo "$line" | jq -r '.genome_assembly')
#   file_name="${download_directory}/ENCODE/ZNF507/${genome_assembly}.$(basename "$file_url")"
#   curl -o "$file_name" -J -L "https://www.encodeproject.org${file_url}"
#   gunzip -c "$file_name" > "${file_name%.bed.gz}.narrowPeak.bed" && rm "$file_name"
# done <<< "$file_info"


#SOX6-K562
#To be stored in data/SOX6_K562 folder
wget https://www.encodeproject.org/files/ENCFF002DPJ/@@download/ENCFF002DPJ.fastq.gz #Replicate 1; Read 1
wget https://www.encodeproject.org/files/ENCFF002EGF/@@download/ENCFF002EGF.fastq.gz #Replicate 1; Read 2
wget https://www.encodeproject.org/files/ENCFF002DUA/@@download/ENCFF002DUA.fastq.gz #Replicate 2; Read 1
wget https://www.encodeproject.org/files/ENCFF002DPM/@@download/ENCFF002DPM.fastq.gz #Replicate 2; Read 2
#Download control ENCSR173USI
wget https://www.encodeproject.org/files/ENCFF002EFF/@@download/ENCFF002EFF.fastq.gz #Replicate 1; Read 1
wget https://www.encodeproject.org/files/ENCFF002EFH/@@download/ENCFF002EFH.fastq.gz #Replicate 1; Read 2
wget https://www.encodeproject.org/files/ENCFF002EFD/@@download/ENCFF002EFD.fastq.gz #Replicate 2; Read 1
wget https://www.encodeproject.org/files/ENCFF002EFA/@@download/ENCFF002EFA.fastq.gz #Replicate 2; Read 2
#Download ENCODE processed data file
folder=ENCODE/SOX6
mkdir -p "$download_directory/$folder" #in case the folder is not created yet
download_ENCODE ENCSR788RSW $folder

#ZNF23-eGFP-HEK293
#To be stored in data/ZNF23_eGPT_HEK293 folder
wget https://www.encodeproject.org/files/ENCFF311MID/@@download/ENCFF311MID.fastq.gz #Replicate 1; Read 1
wget https://www.encodeproject.org/files/ENCFF258FCV/@@download/ENCFF258FCV.fastq.gz #Replicate 1; Read 2
wget https://www.encodeproject.org/files/ENCFF767STU/@@download/ENCFF767STU.fastq.gz #Replicate 2; Read 1
wget https://www.encodeproject.org/files/ENCFF783NMT/@@download/ENCFF783NMT.fastq.gz #Replicate 2; Read 2
#Download control ENCSR007GUS
wget https://www.encodeproject.org/files/ENCFF873ZZU/@@download/ENCFF873ZZU.fastq.gz #Read 1
wget https://www.encodeproject.org/files/ENCFF611URR/@@download/ENCFF611URR.fastq.gz #Read 2
#Download ENCODE processed data file
folder=ENCODE/ZNF23
mkdir -p "$download_directory/$folder" #in case the folder is not created yet
download_ENCODE ENCSR679KXF $folder

###########################################################################################
########### FOR META-ANALYSIS #############################################################
# ENCODE
gene=146 #ZNF146
meta_directory=$wd/meta_analysis_files/$gene
mkdir -p $meta_directory
# hg19
# Unique: Fold change over control BigWig files
wget -O $meta_directory/ENCFF448SWT.hg19.bigWig https://www.encodeproject.org/files/ENCFF448SWT/@@download/ENCFF448SWT.bigWig
# Primary alignment files (BAM) - unfiltered
wget -O $meta_directory/ENCFF912YSO_rep1.hg19.bam https://www.encodeproject.org/files/ENCFF912YSO/@@download/ENCFF912YSO.bam
wget -O $meta_directory/ENCFF898DVS_rep2.hg19.bam https://www.encodeproject.org/files/ENCFF898DVS/@@download/ENCFF898DVS.bam
# hg38
# Unique: Fold change over control BigWig files
wget -O $meta_directory/ENCFF887JZO.hg38.bigWig https://www.encodeproject.org/files/ENCFF887JZO/@@download/ENCFF887JZO.bigWig
# Primary alignment files (BAM) - unfiltered
wget -O $meta_directory/ENCFF652WXQ_rep1.hg38.bam https://www.encodeproject.org/files/ENCFF652WXQ/@@download/ENCFF652WXQ.bam
wget -O $meta_directory/ENCFF830UTZ_rep2.hg38.bam https://www.encodeproject.org/files/ENCFF830UTZ/@@download/ENCFF830UTZ.bam

gene=507 #ZNF507
meta_directory=$wd/meta_analysis_files/$gene
mkdir -p $meta_directory
# hg19
# Unique: Fold change over control BigWig files
wget -O $meta_directory/ENCFF441LZN.hg19.bigWig https://www.encodeproject.org/files/ENCFF441LZN/@@download/ENCFF441LZN.bigWig
# Primary alignment files (BAM) - unfiltered
wget -O $meta_directory/ENCFF579IIE_rep1.hg19.bam https://www.encodeproject.org/files/ENCFF579IIE/@@download/ENCFF579IIE.bam
wget -O $meta_directory/ENCFF995EGD_rep2.hg19.bam https://www.encodeproject.org/files/ENCFF995EGD/@@download/ENCFF995EGD.bam
#hg38
# Unique: Fold change over control BigWig files
wget -O $meta_directory/ENCFF021QXY.hg38.bigWig https://www.encodeproject.org/files/ENCFF021QXY/@@download/ENCFF021QXY.bigWig
# Primary alignment files (BAM) - unfiltered
wget -O $meta_directory/ENCFF156HSL_rep1.hg38.bam https://www.encodeproject.org/files/ENCFF156HSL/@@download/ENCFF156HSL.bam
wget -O $meta_directory/ENCFF582ETC_rep2.hg38.bam https://www.encodeproject.org/files/ENCFF582ETC/@@download/ENCFF582ETC.bam


# Other files needed: 
# L1Base annotation: 
# hg19
wget -P meta_analysis_files -O intact.hg19.bed https://l1base.charite.de/BED/hsflil1_3836.bed
wget -P meta_analysis_files -O disrupted.hg19.bed https://l1base.charite.de/BED/hsorf2l1_3836.bed
# hg38
wget -P meta_analysis_files -O intact.hg38.bed https://l1base.charite.de/BED/hsflil1_8438.bed
wget -P meta_analysis_files -O disrupted.hg38.bed https://l1base.charite.de/BED/hsorf2l1_8438.bed
# RepeatMasker annotation: Genome/hg19.L1.bed
# Mappability track:
# ZNF146-507-Analysis-on-Pangenome/wgEncodeCrgMapabilityAlign100mer.hg19.bigWig #hg19
wget -P -O k100.Umap.MultiTrackMappability.hg38.bw https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k100.Umap.MultiTrackMappability.bw #hg38