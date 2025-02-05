#!/bin/bash
# To load mugqic module
# export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
# module use $MUGQIC_INSTALL_HOME/modulefiles
export wd=/home/hoangnhi/projects/def-bourqueg/hoangnhi/ZNF146-507-Analysis-on-Pangenome
export results_dir=$wd/enrichment_calc # all results from this analysis will be put in this dir
mkdir -p $results_dir
cd $wd
module load StdEnv/2020
module load bedtools/2.30.0

# Need to unload mugqic's python module because of errors when running IDR
module_name=mugqic/python/3.10.4
if module list 2>&1 | grep -q "$module_name"; then
    module unload "$module_name"
fi


# Generate a list of directories of result files created from GenPipes
##################################################################################
############################# ENRICHMENT CALCULATION #############################
enrichment_calc() {
    local bed_path=$1   # to locate the bed file
    local bed_file=$2   # to take the file name and create downstream files
    local csv_file=$3   # to write the enrichment csv file
    local ref=$4        # to use which ref genome
    local blacklist=$5  # to use which blacklist file
    local row_title=$6  # the title of each calculation (1st column)

    ref_l1="$wd/Genome/$ref.L1.bed" # repeatmasker file
    rm=$(wc -l $ref_l1) # total number of L1 peaks in the repeatmasker file
    file=$bed_path/$bed_file
    echo $file
    peaks=$(wc -l "$file.narrowPeak.bed")
    sort -k1,1 -k2,2n "$file.narrowPeak.bed" > "$file.sorted.bed"
    if ! grep -Eq '^(chr)' "$file.sorted.bed"; then
        awk '{print "chr"$0}' "$file.sorted.bed" > temp_file && mv temp_file "$file.sorted.bed"
        # awk -i inplace '{print "chr"$0}' "$file.sorted.bed" # apply only to hg19 GenPipes file because somehow they don't includr "chr"
    fi

    # count overlap with L1 and calculate p_obs
    # NEED TO OPTIMIZE: I should've count the peaks only instead of creating a new bed file
    bedtools intersect -a "$file.sorted.bed" -b $ref_l1 -u > "$bed_path/all_intersection.$ref.bed"
    obs=$(wc -l "$bed_path/all_intersection.$ref.bed") # total number of peaks overlapping with L1
    p_obs=$(echo ${obs%% *} / ${rm%% *} | bc -l) # proportion of peaks overlapping with L1

    p_shuffle=0 
    # shuffle 10 times, exclude blacklist region (not so significant but the paper method said so)
    for _ in {1..10}
    do 
        bedtools shuffle -i "$file.sorted.bed" -g "Genome/human.$ref._noCHR.genome" -noOverlapping -maxTries 1000 -excl $blacklist > "$file.shuffle.bed"
        awk '{print "chr"$0}' "$file.shuffle.bed" > temp_file && mv temp_file "$file.shuffle.bed"
        # awk -i inplace '{print "chr"$0}' "$file.shuffle.bed" # output from bedtools doesn't have "chr"$0
        sort -k1,1 -k2,2n "$file.shuffle.bed" -o "$file.shuffle.bed"
        bedtools intersect -a "$file.shuffle.bed" -b $ref_l1 -u > "$bed_path/all_intersection_shuffle.$ref.bed"
        shuffle=$(wc -l "$bed_path/all_intersection_shuffle.$ref.bed") # total number of peaks overlapping with L1 after shuffling
        echo "${shuffle%% *}"
        p_shuffle=$(echo $p_shuffle + ${shuffle%% *} / ${rm%% *} | bc -l) # summary of proportion of peaks overlapping with L1 after shuffling
        echo "$p_shuffle"
    done


    p_shuffle=$(echo $p_shuffle/10 | bc -l) # average of proportion of peaks overlapping with L1 after shuffling

    enr=$(echo $p_obs / $p_shuffle | bc -l)
    echo "enrichment $enr"
    
    echo -e "${row_title}_${gene}_${ref}\t${peaks%% *}\t${obs%% *}\t$enr" >> "$csv_file"
}
# cat 146_vgmap_results/chm13/callpeaks/*_linear_peaks.bed > 146_vgmap_results/chm13/callpeaks/all_peaks.narrowPeak.bed 
bed_path=146_vgmap_results/chm13/callpeaks   # to locate the bed file
bed_file=all_peaks   # to take the file name and create downstream files
csv_file=146_vgmap_results/chm13/callpeaks/enrichment.csv   # to write the enrichment csv file
ref=t2t        # to use which ref genome
blacklist=Genome/Blacklist/t2t.excluderanges.bed  # to use which blacklist file
row_title=146_vgmap_chm13
enrichment_calc $bed_path $bed_file $csv_file $ref $blacklist $row_title


IDR_creator(){
    # Run IDR: Create IDR directory
    mkdir -p "$path_idr/IDR"
    cd "$path_idr/IDR" || return
    
    # Create a pooled-replicate narrowPeak file like ENCODE (*IMPORTANT: it's the .narrowPeak file, not the .bed file)
    narrowPeak1="${directory[0]/#./$wd}/${file_name[0]}.narrowPeak"
    narrowPeak2="${directory[1]/#./$wd}/${file_name[1]}.narrowPeak"
    cat $narrowPeak1 $narrowPeak2 > "pooled_file.narrowPeak" #double check here

    # Get peaks passing IDR threshold of 5%
    idr --sample $narrowPeak1 $narrowPeak2 --peak-list "pooled_file.narrowPeak" --input-file-type narrowPeak --output-file "idrValues.txt" --plot
        # IDR output file idrValues.txt
        # Columns 1-10 are same as pooled common peaks narrowPeak columns 
        # Col 11: -log10(local IDR value)
        # Col 12: -log10(global IDR value)    <-------
        # Col 15: ranking measure from Rep1
        # Col 19: ranking measure from Rep2
    echo 'IDR file was created'

    # Calculate the negative log base 10 of 0.05 -log(p)/log(10) to filter peaks in col 12 of idrValues.txt
    IDR_THRESH=0.05
    IDR_THRESH_TRANSFORMED=$(awk -v p=${IDR_THRESH} 'BEGIN{print -log(p)/log(10)}') 
    awk 'BEGIN{OFS="\t"} $12>='"${IDR_THRESH_TRANSFORMED}"' {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' "idrValues.txt" \
        | sort | uniq | sort -k7n,7n > $gene.idr.narrowPeak.bed

    # Number of peaks passing IDR threshold of 5%
    NPEAKS_IDR=$(cat $gene.idr.narrowPeak.bed | wc -l)
    echo "Number of peaks passing IDR threshold of 5% $NPEAKS_IDR"

    # Filter using blacklist
    if [[ "$ref" == "hg19" ]]; then
        awk -i inplace '{print "chr"$0}' "$gene.idr.narrowPeak.bed"
    fi
    bedtools intersect -v -a $gene.idr.narrowPeak.bed -b "$wd/${blacklist_file}" | grep -P 'chr[\dXY]+[ \t]' \
        | awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' > $gene.idr.filt.narrowPeak.bed

    cd $wd || return # go back to working directory
}

# This is results before 18/03 having 2 replications
write_to_csv_2rep(){
    local ref=$1 
    local gene=$2
    local blacklist_file=$3
    local dir1=$4
    local dir2=$5
    local csv_file=$6
    local IDR_check=$7 # 0 is to calculate IDR, 1 is to skip and use the normal IDR file, 2 is to skip and use the 95% IDR file

    
    directory=("$(dirname $dir1)" "$(dirname $dir2)")
    local path_idr=$(echo "$dir1" | awk -F'/' '{OFS="/"; for (i=1; i<=3; i++) printf "%s%s", $i, (i<3 ? OFS : ""); print ""}')
    
    file_name=("$(basename "${dir1%.narrowPeak.bed}")" "$(basename "${dir2%.narrowPeak.bed}")")
    ref=$(echo "$dir1" | awk -F'/' '{print $2}' | cut -d'_' -f3)
    gene=$(echo "$dir1" | awk -F'/' '{print $2}' | cut -d'_' -f1)
    blacklist_file=$(find Genome/Blacklist -type f -name "*${ref#hg}*")

    # Enrichment of ENCODE file if reference genome is hg19 or hg38
    if [[ "$ref" == "hg19" || "$ref" == "hg38" ]]; then
        encode_directory="$wd/GenPipes_sets/data/ENCODE/$gene"
        encode_file=$(basename "$(find $encode_directory -type f -name "*${ref#hg}*.narrowPeak.bed")")
        enrichment_calc "$encode_directory" "${encode_file%.narrowPeak.bed}" $csv_file $ref $blacklist_file "ENCODE"
    fi

    # Enrichment of each of 2 replications
    # enrichment_calc "${directory[0]}" "${file_name[0]}" $csv_file $ref $blacklist_file "$gene/Rep1"
    # enrichment_calc "${directory[1]}" "${file_name[1]}" $csv_file $ref $blacklist_file "$gene/Rep2"

    # IDR_creator
    if [ "$IDR_check" == 0 ]; then
        IDR_creator
    fi

    # Enrichment of IDR
    enrichment_calc "$path_idr/IDR" "$gene.idr.filt" $csv_file $ref $blacklist_file "GenPipes_IDR"

    if [ "$IDR_check" == 2  ]; then
        enrichment_calc "$path_idr/IDR" "$gene.90percent.idr.filt" $csv_file $ref $blacklist_file "GenPipes_IDR_95%"
    fi

    echo "" >> "$csv_file"
}

#### FORCE ALL PEAKS TO HAVE THE SAME SIZE 500BP ####
write_to_csv_500bp(){
    local ref=$1 
    local gene=$2
    local blacklist_file=$3
    local dir1=$4
    local dir2=$5
    local csv_file=$6

    # Create 500bp peaks bed file for IDR (ENCODE)
    encode_500bp_directory=$wd/GenPipes_sets/data/ENCODE_500bp/$gene
    mkdir -p "$encode_500bp_directory"
    if [[ "$ref" == "hg19" || "$ref" == "hg38" ]]; then
        encode_directory="$wd/GenPipes_sets/data/ENCODE/$gene"
        encode_file=$(basename "$(find $encode_directory -type f -name "*${ref#hg}*.narrowPeak.bed")")
        echo "$encode_file"
        awk 'BEGIN {OFS="\t"} {mid = int(($3 + $2) / 2); start = mid - 250; end = mid + 250; print $1, start, end, $4, $5, $6, $7, $8, $9, $10}'\
            "$encode_directory/$encode_file" > "$encode_500bp_directory/${ref#hg}.$gene.500bp.narrowPeak.bed"
        enrichment_calc "$encode_500bp_directory" "${ref#hg}.$gene.500bp" $csv_file $ref $blacklist_file "ENCODE 500bp"
    fi

    # Create 500bp peaks bed file for IDR (GenPipes)
    path_idr=$(echo "$dir1" | awk -F'/' '{OFS="/"; for (i=1; i<=3; i++) printf "%s%s", $i, (i<3 ? OFS : ""); print ""}')/IDR
    echo $path_idr
    awk 'BEGIN {OFS="\t"} {mid = int(($3 + $2) / 2); start = mid - 250; end = mid + 250; print $1, start, end, $4, $5, $6, $7, $8, $9, $10}'\
        "$path_idr/$gene.idr.filt.narrowPeak.bed" > "$path_idr/$gene.idr.filt.500bp.narrowPeak.bed"
    enrichment_calc "$path_idr" "$gene.idr.filt.500bp" $csv_file $ref $blacklist_file "GenPipes 500bp"
    echo "" >> "$csv_file"
}

# This is results after 18/03 having 1 replication - i.e. "_design2_" in the name
write_to_csv_1rep(){
    local ref=$1 
    local gene=$2
    local blacklist_file=$3
    local dir1=$4
    local csv_file=$5
    
    file_name=("$(basename "${dir1%.narrowPeak.bed}")")

    # Enrichment of ENCODE file if reference genome is hg19 or hg38
    if [[ "$ref" == "hg19" || "$ref" == "hg38" ]]; then
        encode_directory="$wd/GenPipes_sets/data/ENCODE/$gene"
        encode_file=$(basename "$(find $encode_directory -type f -name "*${ref#hg}*.narrowPeak.bed")")
        enrichment_calc "$encode_directory" "${encode_file%.narrowPeak.bed}" $csv_file $ref $blacklist_file "ENCODE"
    fi

    # Enrichment of GenPipes file
    enrichment_calc $(dirname $dir1) $file_name $csv_file $ref $blacklist_file "GenPipes"
}

type="design2"
ref="t2t"
# directory_list=($(find $wd -type f -name '*peaks.narrowPeak.bed' | grep -vE '/(report|146_results_24controls_hg19)/'  | sort))
directory_list=($(find $wd -type f -name '*peaks.narrowPeak.bed' | grep -vE '/report/' | grep "_${type}_${ref}"  | sort))
echo "${directory_list[@]}"

# Put all enrichment result files in 1 folder
date=$(date "+%Y-%m-%d")
enr_file_name="${date}_enrichment_calc"
csv_file="$results_dir/$enr_file_name.csv"
echo -e "Experiment\tNumPeaks\tL1Peaks\tEnrichment" > "$csv_file"

########### Analyze only unique T2T regions
t2t_unique_bed=$wd/Genome/hgUnique.hg38.bed

########### Loop for combined replication from GenPipes after 18/3
while IFS= read -r dir1; do
    ref=$(echo "${dir1#$wd}" | awk -F'/' '{print $2}' | cut -d'_' -f3)
    gene=$(echo "${dir1#$wd}" | awk -F'/' '{print $2}' | cut -d'_' -f1)
    blacklist_file=$(find Genome/Blacklist -type f -name "*${ref#hg}*")
    echo $dir1 $ref $gene $blacklist_file
    unique_bed=$(dirname $dir1)/$(basename $dir1 .narrowPeak.bed).t2t.unique.narrowPeak.bed
    bedtools intersect -a $dir1 -b $t2t_unique_bed -u > $unique_bed
    write_to_csv_1rep $ref $gene $blacklist_file $unique_bed $csv_file
done < <(printf "%s\n" "${directory_list[@]}")
dir1="${directory_list[0]}"

########### Loop for 2 replications
# while IFS= read -r dir1 && IFS= read -r dir2; do
#     echo "Directory 1: $dir1"
#     echo "Directory 2: $dir2"

#     ref=$(echo "$dir1" | awk -F'/' '{print $2}' | cut -d'_' -f3)
#     gene=$(echo "$dir1" | awk -F'/' '{print $2}' | cut -d'_' -f1)
#     blacklist_file=$(find Genome/Blacklist -type f -name "*${ref#hg}*")
    
#     #### FORCE ALL PEAKS TO HAVE THE SAME SIZE 500BP ####
#     # write_to_csv_500bp $ref $gene $blacklist_file $dir1 $dir2 $csv_file
#     ###################

#     # write_to_csv $ref $gene $blacklist_file $dir1 $dir2 $csv_file 0
    
# done < <(printf "%s\n" "${directory_list[@]}")

awk 'BEGIN {FS=","; OFS="\t"} {for (i=1; i<=NF; i++) printf "%s%s", $i, (i<NF?OFS:ORS)}' $csv_file > "$results_dir/$enr_file_name.tsv"

####################################################################################################################
############################# TROUBLESHOOTING: OVERLAPPING ENCODE AND GENPIPES RESUTLS #############################
# All the bedtools intersect between ENCODE and GENPIPES will be in enrichment_calc folder
########### A different version of it
directory_list=()
type="design2"
refs=("hg19" "hg38" "t2t")
genes=("146" "507")
for ref in "${refs[@]}"; do
    for gene in "${genes[@]}"; do
        directory_list+=("${gene}_${type}_${ref}")
    done
done
encode_genpipes_dir="$wd/enrichment_calc/ENCODEvsGENPIPES"
mkdir -p $encode_genpipes_dir

add_peaks_to_list(){
    local file=$1
    count=$(wc -l "$file")
    results+=("${count%% *}")
}
overlap(){
    local encode_IDR=$1
    local genpipes_IDR=$2
    local type=$3
    bedtools intersect -a $encode_IDR -b $genpipes_IDR -u > "$encode_genpipes_dir/overlap.$type.$gene.bed"
    bedtools intersect -a $genpipes_IDR -b $encode_IDR -u > "$encode_genpipes_dir/overlap2.$type.$gene.bed"
    bedtools intersect -a $encode_IDR -b $genpipes_IDR -v > "$encode_genpipes_dir/encode_notingenpipes.$type.$gene.bed"
    bedtools intersect -a $genpipes_IDR -b $encode_IDR -v > "$encode_genpipes_dir/genpipes_notinencode.$type.$gene.bed"
    
    add_peaks_to_list "$encode_genpipes_dir/overlap.$type.$gene.bed"
    add_peaks_to_list "$encode_genpipes_dir/overlap2.$type.$gene.bed"
    add_peaks_to_list "$encode_genpipes_dir/encode_notingenpipes.$type.$gene.bed" 
    add_peaks_to_list "$encode_genpipes_dir/genpipes_notinencode.$type.$gene.bed" 

    results+=(".")
}

# Put all enrichment result files in 1 folder
date=$(date "+%Y-%m-%d")
csv_file="$wd/enrichment_calc/${date}_ENCODEvsGENPIPES.csv"
truncate -s 0 $csv_file
# Write the labels to the csv file first, as a column
my_list=()
for _ in {1..15}; do
    my_list+=(".")
done
label_list=("ENCODE/GENPIPES" "e/g all peaks" "g/e all peaks" "e+/p- all peaks" "e-/p+ all peaks" "." \
            "e/g L1 peaks" "g/e L1 peaks" "e+/p- L1 peaks" "e-/p+ L1 peaks" "." \
            "E peak size" "G peak size") # "." "+/+ Top 15")
#label_list=("${label_list[@]}" "${my_list[@]}" "+/- Top 15" "${my_list[@]}" "-/+ Top 15" "${my_list[@]}")
for item in "${label_list[@]}"; do
    echo "$item" >> "$csv_file"
done   

# ok=true
while IFS= read -r dir1 && IFS= read -r dir2; do
    echo "Directory 1: $dir1"
    echo "Directory 2: $dir2"
    # directory=("$(dirname $dir1)" "$(dirname $dir2)")

    for dir in "$dir1" "$dir2"; do
        path_idr="$wd/$dir/peak_call/IDR"

        ref=$(echo "$dir" | cut -d'_' -f3)
        gene=$(echo "$dir"| cut -d'_' -f1)

        results=("$gene.$ref")

        # Overlap 2 narrowpeak files
        encode_directory="$wd/GenPipes_sets/data/ENCODE/$gene"
        encode_all=$(basename "$(find $encode_directory -type f -name "*${ref#hg}*.sorted.bed")")
        genpipes_all="$path_idr/$gene.idr.filt.sorted.bed"

        overlap "$encode_directory/$encode_all" "$genpipes_all" "all" #results



        # Overlap 2 L1 files
        encode_IDR="$wd/GenPipes_sets/data/ENCODE/$gene/all_intersection.$ref.bed"
        genpipes_IDR="$path_idr/all_intersection.$ref.bed"
        bedtools intersect -a "$encode_directory/$encode_all" -b "$wd/Genome/$ref.L1.bed" -u > $encode_IDR
        bedtools intersect -a "$genpipes_all" -b "$wd/Genome/$ref.L1.bed" -u > $genpipes_IDR
        
        overlap "$encode_IDR" "$genpipes_IDR" "L1" #results

        results+=("$(awk '{total += $3 - $2; count++} END {print total/count}' "$encode_directory/$encode_all")" \
                  "$(awk '{total += $3 - $2; count++} END {print total/count}' $genpipes_all)")

        # I don't need this step yet, will work on it later!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # List top 15 L1 subclasses in overlap
        # ref_l1="$wd/Genome/$ref.L1.bed"
        # top15_files=("$encode_genpipes_dir/overlap.L1.$gene.bed" "$encode_genpipes_dir/encode_notingenpipes.L1.$gene.bed"\
        #              "$encode_genpipes_dir/genpipes_notinencode.L1.$gene.bed")
        # for item in "${top15_files[@]}"; do
        #     count=0
        #     while IFS2= read -r line; do
        #         results+=("$line")
        #         ((count++))
        #     done < <(bedtools intersect -wb -a $item -b $ref_l1 | cut -f14 | sort | uniq -c | sort -nr | head -n 15)
        #     echo $count
        #     for ((i=0; i<16-$count; i++)); do
        #         results+=(".")
        #     done
        # done

        echo -e "${results[@]}"

        temp_file="$wd/enrichment_calc/temp.csv"
        truncate -s 0 $temp_file

        index=0
        # Read the existing CSV file line by line
        while IFS2= read -r line; do
            # Append the new list element to the current line
            echo -e "${line}\t${results[$index]}" >> "$temp_file"
            # Increment the index
            ((index++))
        done < "$csv_file"
        # Replace the old CSV file with the new temporary file
        mv "$temp_file" "$csv_file"
    # ok=false
    done
done < <(printf "%s\n" "${directory_list[@]}")
awk 'BEGIN {FS=","; OFS="\t"} {for (i=1; i<=NF; i++) printf "%s%s", $i, (i<NF?OFS:ORS)}' $csv_file > "$wd/enrichment_calc/${date}_ENCODEvsGENPIPES.tsv"

#####################################################################################################################
############################# Filter GenPipes too large peaks then calculate enrichment again #############################
directory_list=("146_results_hg19" "146_results_hg38" "507_results_hg19" "507_results_hg38")
upper_bound=("1069" "1063" "1869" "1773") # from R plots
date=$(date "+%Y-%m-%d")
csv_file="$wd/enrichment_calc/${date}_enrichment_calc_90percent.csv"
echo -e "Rep\tNumPeaks\tL1Peaks\tEnrichment" > "$csv_file"

for ((i=0; i<${#directory_list[@]}; i++)); do
    dir=${directory_list[i]}
    ref=$(echo "$dir" | cut -d'_' -f3)
    gene=$(echo "$dir"| cut -d'_' -f1)
    blacklist_file=$(find Genome/Blacklist -type f -name "*${ref#hg}*")
    path_idr="$wd/$dir/peak_call/IDR"
    input_bed_file="$path_idr/$gene.idr.filt.narrowPeak.bed"
    output_bed_file="$path_idr/$gene.90percent.idr.filt.narrowPeak.bed"
    awk -v ub="${upper_bound[i]}" '{ if (($3-$2) < ub) print $0; }' "$input_bed_file" > "$output_bed_file"
    enrichment_calc $path_idr "$gene.90percent.idr.filt" $csv_file $ref $blacklist_file 1
done

#####################################################################################################################
############################# CREATE PIE CHART for ZNF146 and ZNF507 (hg19, 38 and t2t) #############################
cd $wd
# directory_list=($(find . -type f -name '*peaks.narrowPeak.bed' | grep -vE '/(report|146_results_24controls_hg19)/' | sort))
directory_list=($(find $wd -type f -name '*peaks.narrowPeak.bed' | grep -vE '/report/' | grep '_design2_'  | sort))
genes=("146" "507")
hg=("hg19" "hg38")

L1_peaks_count(){
    local gene=$1
    local ref=$2
    local path=$3
    local line1=$(bedtools intersect -a "$path" -b "$wd/Genome/$ref.L1.bed" -u | wc -l)
    local all=$(wc -l "$path")
    echo ${line1%% *} / ${all%% *}*100 | bc -l
}

piechart_calc(){
    local ref=$1
    local path146=$2    # directory to 146 peaksfile
    local path507=$3    # directory to 507 peaks file
    local name=$4       # to write in the 1st column
    local csv_file=$5
    local list=("$name")

    # First, calculate the percentage of L1 peaks for each gene
    list+=("$(L1_peaks_count "146" "$ref" "$path146")" "$(L1_peaks_count "507" "$ref" "$path507")")

    # Second, overlap 2 L1 peaks to see that ZNF146 and 507 don't target the same location
    count_146=$(wc -l "$path146")
    count_507=$(wc -l "$path507")
    inter=$(bedtools intersect -a "$path146" -b "$path507" -u | wc -l)
    noninter_146=$(echo ${count_146%% *} - ${inter%% *} | bc -l)
    noninter_507=$(echo ${count_507%% *} - ${inter%% *} | bc -l)
    list+=("$noninter_146" "${inter%% *}" "$noninter_507")
    echo "${list[@]}" | tr ' ' '\t' >> $csv_file
}

date=$(date "+%Y-%m-%d")
csv_file="$wd/enrichment_calc/${date}_piechart.csv"
echo -e "Reference\t146_L1Peaks\t507_L1Peaks\t146_noninter\tInter\t507_non_inter" > "$csv_file" # 6 columns
#control
for ref in "${hg[@]}"; do
    # Do calculations for ENCODE files
    if [[ "$ref" == "hg19" || "$ref" == "hg38" ]]; then
        encode_directories=($(for gene in "${genes[@]}"; do echo "$wd/GenPipes_sets/data/ENCODE/$gene"; done))
        encode_files=($(for edir in "${encode_directories[@]}"; do echo "$edir/$(basename "$(find "$edir" -type f -name "*${ref#hg}*.sorted.bed")")"; done))
        piechart_calc "$ref" "${encode_files[0]}" "${encode_files[1]}" "$ref.encode" "$csv_file"
    fi

    # Do calculations for GenPipes files
    genpipes_files=("$wd/146_${type}_${ref}/peak_call/146/ZNF146/146.ZNF146_peaks.sorted.bed" "$wd/507_${type}_${ref}/peak_call/507/ZNF507/507.ZNF507_peaks.sorted.bed")
    piechart_calc "$ref" "${genpipes_files[0]}" "${genpipes_files[1]}" "$ref.genpipes" "$csv_file"
done
awk 'BEGIN {FS=","; OFS="\t"} {for (i=1; i<=NF; i++) printf "%s%s", $i, (i<NF?OFS:ORS)}' $csv_file > "$wd/enrichment_calc/${date}_piechart.tsv"

#####################################################################################################################
############################# FORCE ALL PEAKS TO HAVE THE SAME SIZE 500BP #############################

awk 'BEGIN {OFS="\t"} {mid = int(($3 + $2) / 2); start = mid - 250; end = mid + 250; print $1, start, end, $4, $5, $6, $7, $8, $9, $10, $11, $12}' 146_design_hg19/peak_call/IDR/146.idr.filt.sorted.bed > 146_design_hg19/peak_call/IDR/146.idr.filt.sorted.500bp.bed
awk '($3 - $2) > 600' $wd/146_design2_hg19/peak_call/146/ZNF146/all_intersection.hg19.bed > $wd/146_design2_hg19/peak_call/146/ZNF146/all_intersection.hg19.600.bed

bedtools intersect -a /home/hoangnhi/projects/def-bourqueg/hoangnhi/ZNF146-507-Analysis-on-Pangenome/146_design_hg19/peak_call/IDR/146.idr.filt.sorted.1kb.bed -b "$wd/Genome/hg19.L1.bed" -u > "/home/hoangnhi/projects/def-bourqueg/hoangnhi/ZNF146-507-Analysis-on-Pangenome/146_design_hg19/peak_call/IDR/all_intersection.19.1kb.bed"