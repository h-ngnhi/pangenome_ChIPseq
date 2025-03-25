#!/bin/bash 
export wd=$(pwd)
vg_giraffe_graph() {            # To create graph and index for vg giraffe
    local variant_vcf=$1
    local ref=$2
    local graph_type=$3         # Name of the graph
    graph_dir=$wd/Pangenomes/vg_giraffe/$graph_type
    mkdir -p $graph_dir
    cd $graph_dir || exit
    if [ "$graph_type" == "chm13" ]; then
        vg autoindex --workflow giraffe -r $ref -p $graph_type --threads 64
    else
        vg autoindex --workflow giraffe -r $ref -v $variant_vcf -p $graph_type --threads 64
    fi
    mv $graph_type.giraffe.gbz $graph_type.gbz
    cd $wd || exit
}


vg_map_graph() {                 # To create graph and index for vg map
    local variant_vcf=$1
    local ref=$2
    local graph_type=$3         # Name of the graph

    graph_dir=$wd/Pangenomes/vg_map/$graph_type
    echo $graph_dir
    mkdir -p $graph_dir
    mkdir -p $graph_dir/graphs
    cd $graph_dir || exit
    if [ "$graph_type" == "chm13" ]; then
        (seq 1 22; echo X; echo M) | parallel -j 6 "vg construct -S -a -p -C -R chr{} -r $ref -t 1 -m 32 > graphs/chr{}.vg"
        vg ids -j -m mapping.txt graphs/chr*.vg
        vg index -x graph.xg graphs/chr*.vg -t 40
        vg index -g graph.gcsa -f mapping.txt graphs/chr*.vg -t 40
    else
        module load tabix
        mkdir -p vcf
        (seq 1 22; echo X; echo M) | parallel -j 24 "tabix -f -h $variant_vcf chr{} > vcf/chr{}.vcf ; bgzip vcf/chr{}.vcf ; tabix vcf/chr{}.vcf.gz"
        (seq 1 22; echo X; echo M) | parallel -j 6 "vg construct -S -a -p -C -R chr{} -v vcf/chr{}.vcf.gz -r $ref -t 1 -m 32 > graphs/chr{}.vg"
        vg ids -j -m mapping.txt graphs/chr*.vg
        # Adding variants --> need to create gbwt (haplotype-aware index) file
        (seq 1 22; echo X; echo M) | parallel -j 8 "vg gbwt -x graphs/chr{}.vg -v vcf/chr{}.vcf.gz -o graphs/chr{}.gbwt"
        vg gbwt -m -o graph.gbwt graphs/chr*.gbwt
        vg index -x graph.xg graphs/chr*.vg 
        # Prune the graph for gcas index
        mkdir -p graphs_pruned
        (seq 1 22; echo X; echo M) | parallel -j 6 "vg prune -u -a -m mapping.txt -g graphs/chr{}.gbwt graphs/chr{}.vg > graphs_pruned/chr{}.pruned.vg"
        vg index -g graph.gcsa -f mapping.txt -Z 4096 graphs_pruned/*.vg -t 40
    fi
    cd $wd || exit
}

vg_map_convert() {                  # To index from gbz graph for vg_map
    local ref=$1
    local graph_dir=$wd/Pangenomes/vg_map/${ref}_converted
    local gbz_dir=$wd/Pangenomes/vg_giraffe/$ref

    mkdir -p $graph_dir
    cd $graph_dir

    # Adding variants --> need to create gbwt (haplotype-aware index) file
    vg gbwt -Z $gbz_dir/$ref.gbz -o graph.gbwt
    vg convert $gbz_dir/$ref.gbz > graph.vg
    vg ids -j -m mapping.txt graph.vg
    vg prune -u -a -m mapping.txt -g graph.gbwt graph.vg > graph.pruned.vg
    vg index -x graph.xg graph.pruned.vg
    vg index -g graph.gcsa -f mapping.txt -Z 4096 graph.pruned.vg -t 40
}

graph_bam() {
    local results_dir=$1

    cd $results_dir || exit
    module load samtools
    samtools view treatment_alignments.bam | awk '{print $5}' > mapq_scores.txt
    awk '$1 ==0 {count++} END {print "MAPQ = 0:", count+0}' mapq_scores.txt
    awk '$1 >0 && $1 < 30 {count++} END {print "0 < MAPQ < 30:", count+0}' mapq_scores.txt
    awk '$1 >= 30 && $1 < 60 {count++} END {print "30 <= MAPQ < 60:", count+0}' mapq_scores.txt
    awk '$1 == 60 {count++} END {print "MAPQ = 60:", count+0}' mapq_scores.txt
    cd $wd || exit
}

alignment() {                   # To map reads using vgmap
    local pipeline=$1       # vg_giraffe or vg_map
    local graph_type=$2     # chm13 or vcfbub
    local results_dir=$3    # result directory
    local forward=($4)      # Array of forward reads
    local reverse=($5)      # Array of reverse reads
    local input=$6          # treatment or control
    local graph_dir=$wd/Pangenomes/$pipeline/$graph_type
    mkdir -p $results_dir
    cd $results_dir || exit
    alignments=()
    if [ "$pipeline" == "vg_giraffe" ]; then
        for i in "${!forward[@]}"; do
            vg giraffe -Z $graph_dir/$graph_type.gbz -m $graph_dir/$graph_type.min -d $graph_dir/$graph_type.dist -f ${forward[i]} -f ${reverse[i]} -t 64 > ${input}_alignments_$((i+1)).gam
            alignments+=("${input}_alignments_$((i+1)).gam")
        done
        cat "${alignments[@]}" > ${input}_alignments.gam
        vg surject -x $graph_dir/$graph_type.gbz -b ${input}_alignments.gam > ${input}_alignments.bam
    elif [ "$pipeline" == "vg_map" ]; then
        for i in "${!forward[@]}"; do
            vg map -x $graph_dir/graph.xg -g $graph_dir/graph.gcsa -f ${forward[i]} -f ${reverse[i]} -t 40 -u 1 -m 1 > ${input}_alignments_$((i+1)).gam
            alignments+=("${input}_alignments_$((i+1)).gam")
        done
        cat "${alignments[@]}" > ${input}_alignments.gam
        vg surject -x $graph_dir/graph.xg -b ${input}_alignments.gam > ${input}_alignments.bam
    fi
    vg filter ${input}_alignments.gam -P -fu -q 30 -t 64 > ${input}_alignments.filtered.gam
    vg view -aj ${input}_alignments.filtered.gam > ${input}_alignments.filtered.json
    if [ "$input" == "treatment" ]; then
        graph_bam $results_dir
    fi
    cd $wd || exit
}


split_graph() {                         # To split graphs by chromosomes for GP
    local pipeline=$1                   # vg_giraffe or vg_map
    local graph=$2                      # chm13 or vcfbub or hprc-v1.1-mc-chm13 or L1_vcfbub
    local gf_sif=$3                     # directory to gf image

    graph_dir=$wd/Pangenomes/$pipeline/$graph
    if [ "$pipeline" == "vg_giraffe" ]; then
        #Split the graph into chromosomes
        mkdir -p $graph_dir/graphs
        vg chunk -C -x $graph_dir/$graph.gbz -b $graph_dir/graphs/

        # Change the name of the graph files to the chromosome name
        for file in $graph_dir/graphs/*.vg; do
            # chr_name=$(vg paths -L -x $file | grep 'CHM13#0#chr' | head -n 1 | awk -F'#' '{print $3}' | awk -F'[' '{print $1}')
            chr_name=$(vg paths -L -x $file | head -n 1)
            mv $file $graph_dir/graphs/${chr_name}.vg
        done
    fi

    graph_dir=$graph_dir/graphs
    # Convert graphs to json
    (seq 1 22; echo X; echo Y; echo M) | parallel -j 2 "vg view -j $graph_dir/chr{}.vg > $graph_dir/chr{}.json ; vg stats -r $graph_dir/chr{}.vg | cut -f 2 > $graph_dir/node_range_chr{}.txt"

    # Create ob_graph files
    SIF_IMAGE="$gf_sif"
    # BIND_DIR="/lustre06/project/6002326/hoangnhi/ZNF146-507-Analysis-on-Pangenome/${graph_dir#$wd/}:/mnt_data"
    BIND_DIR="$graph_dir:/mnt_data"

    module load apptainer
    apptainer exec --contain --cleanenv --bind $BIND_DIR $SIF_IMAGE /bin/bash -c "
        chromosomes=\"chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrM\"
        chromosomes=\$(echo \$chromosomes | tr ',' ' ')
        for chromosome in \$chromosomes; do
            graph_peak_caller create_ob_graph /mnt_data/\$chromosome.json
            graph_peak_caller find_linear_path -g /mnt_data/\$chromosome.nobg /mnt_data/\$chromosome.json \$chromosome /mnt_data/\${chromosome}_linear_path.interval
        done
    
    "
}

callpeaks() {
    local treatment=$1              # treatment json file
    local control=$2                # control json file (could be empty if no control)
    local results_dir=$3            # result directory
    local graph=$4                  # chm13 or vcfbub
    local pipeline=$5               # vg_giraffe or vg_map
    local steps=$6                  # steps to run
    local server=$7                 # narval or ashbi
    local gp_sif=$8

    module load apptainer
    if [[ " $steps " =~ 4 ]]; then
        fragment_length() {
            grep -i "predicted fragment length is" | head -n1 | sed -E 's/.*predicted fragment length is[[:space:]]*([0-9]+).*/\1/'
        }
        if [ "$server" == "narval" ]; then
            SIF_IMAGE="$wd/tools/macs2_latest.sif"
            BIND_DIR="$results_dir:/mnt"
            frag_len=$(apptainer exec --contain --cleanenv --bind $BIND_DIR $SIF_IMAGE macs2 predictd -i /mnt/treatment_alignments.bam 2>&1 | fragment_length)
        elif [ "$server" == "ashbi" ]; then
            module load macs2
            frag_len=$(macs2 predictd -i $results_dir/treatment_alignments.bam 2>&1 | fragment_length)
        fi
        read_len=$(zcat $forward_trm | head -2 | tail -1 | wc -c)
        unique_reads=$(grep -Po '"sequence": "\K([ACGTNacgtn]{20,})"' $results_dir/treatment_alignments.filtered.json | sort | uniq | wc -l)

        if (( frag_len < read_len )); then
            # Round up read_len to the nearest multiple of 5.
            frag_len_org=$frag_len
            frag_len=$(( ((read_len + 4) / 5) * 5 ))
        fi
        
    # Prepare the parameters file for using inside the container (no indent for this command)
cat <<EOL > $results_dir/params.txt
fragment_length=$frag_len_org
fragment_length=$frag_len
read_length=$read_len
unique_reads=$unique_reads
graph_dir=Pangenomes/$pipeline/$graph/graphs
treatment=$treatment
control=$control
steps="$steps"
EOL
    else 
        head -n -1 $results_dir/params.txt > tmp && echo steps=\""$steps"\" >> tmp && mv tmp $results_dir/params.txt
    fi

    SIF_IMAGE="$gp_sif"
    # BIND_DIR="/lustre06/project/6002326/hoangnhi/ZNF146-507-Analysis-on-Pangenome:/mnt_data"
    BIND_DIR="$wd:/mnt"
    BIND_DIR_2="$results_dir:/mnt_results"

    apptainer exec --contain --cleanenv --bind $BIND_DIR --bind $BIND_DIR_2 $SIF_IMAGE /bin/bash -c "
        echo_slurm() {
            local message="\$1"
            echo \$message
            echo \$message >&2
        }
        
        while IFS= read -r line; do 
            eval \"export \$line\"
        done < /mnt_results/params.txt
        results_dir=/mnt_results
        graph_dir=/mnt/\$graph_dir

        chromosomes=\"chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrM\"
            
        # Split the treatment and control json files into chromosomes
        # at this step, it is required to pass in all the chromosomes at ones to the command
        if [[ \" \$steps \" =~ 5 ]]; then
            echo_slurm \"Split the treatment and control json files into chromosomes\"
            graph_peak_caller split_vg_json_reads_into_chromosomes \$chromosomes \$results_dir/\$treatment \$graph_dir/
            if [ -n \"\$control\" ]; then
                graph_peak_caller split_vg_json_reads_into_chromosomes \$chromosomes \$results_dir/\$control \$graph_dir/
            fi
            mkdir -p \$results_dir/reads_by_chrom
            mv \$results_dir/*chr*.json \$results_dir/reads_by_chrom/
        fi

        # Call peaks
        if [[ \" \$steps \" =~ 6 ]]; then
            echo_slurm \"Call peaks\"
            mkdir -p \$results_dir/callpeaks/
            treatment_chr=\$(basename \${treatment%.*})
            echo_slurm \$treatment_chr
            ctl_chr=\$(basename \${control%.*})
            echo_slurm \$ctl_chr

            chromosomes=\$(echo \$chromosomes | tr ',' ' ')
            for chromosome in \$chromosomes; do
                echo_slurm \$chromosome
                if [ -n \"\$control\" ]; then
                    graph_peak_caller callpeaks -g \$graph_dir/\$chromosome.nobg -s \$results_dir/reads_by_chrom/\${treatment_chr}_\$chromosome.json -c \$results_dir/reads_by_chrom/\${ctl_chr}_\$chromosome.json -f \$fragment_length -r \$read_length -u \$unique_reads -p True -G 3100000000 -n \$results_dir/callpeaks/\${chromosome}_
                else
                    graph_peak_caller callpeaks -g \$graph_dir/\$chromosome.nobg -s \$results_dir/reads_by_chrom/\${treatment_chr}_\$chromosome.json -f \$fragment_length -r \$read_length -u \$unique_reads -p True -G 3100000000 -n \$results_dir/callpeaks/\${chromosome}_
                fi
            done

            cd \$results_dir/callpeaks/
            echo_slurm \"Call peaks for the whole genome\"
            for chromosome in \$chromosomes; do
                echo_slurm \$chromosome
                graph_peak_caller callpeaks_whole_genome_from_p_values -d \$graph_dir/ -n \"\" -f \$fragment_length -r \$read_length \$chromosome
                graph_peak_caller peaks_to_linear \${chromosome}_max_paths.intervalcollection \$graph_dir/\${chromosome}_linear_path.interval \$chromosome \${chromosome}_linear_peaks.bed
            done

            echo_slurm \"Concatenate_sequence_files\"
            chromosomes=\"chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrM\"
            graph_peak_caller concatenate_sequence_files \$chromosomes all_peaks.fasta
            wc -l all_peaks.intervalcollection
        fi

        if [[ \" \$steps \" =~ 7 ]]; then
            cd \$results_dir/callpeaks/
            echo_slurm \"Call peaks for the whole genome\"
            chromosomes=\$(echo \$chromosomes | tr ',' ' ')
            for chromosome in \$chromosomes; do
                echo_slurm \$chromosome
                if ! ls \$graph_dir/\${chromosome}_linear_path.interval > /dev/null 2>&1; then
                    graph_peak_caller find_linear_path -g \$graph_dir/\$chromosome.nobg \$graph_dir/\$chromosome.json \$chromosome \$graph_dir/\${chromosome}_linear_path.interval
                fi
                graph_peak_caller callpeaks_whole_genome_from_p_values -d \$graph_dir/ -n \"\" -f \$fragment_length -r \$read_length \$chromosome
                graph_peak_caller peaks_to_linear \${chromosome}_max_paths.intervalcollection \$graph_dir/\${chromosome}_linear_path.interval \$chromosome \${chromosome}_linear_peaks.bed
            done

            echo_slurm \"Concatenate_sequence_files\"
            chromosomes=\"chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrM\"
            graph_peak_caller concatenate_sequence_files \$chromosomes all_peaks.fasta
            wc -l all_peaks.intervalcollection
        fi
    "
}

# Editing MAPQ for giraffe distribution
edit_mapq() {
    local markname=$1
    local pipeline=$2
    local ref=$3

    local results_dir=$wd/results/$markname/${pipeline}_${ref}
    awk '
    BEGIN { FS=OFS="\"mapping_quality\":" }
    {
    if (NF > 1 && $2 ~ /^[ \t]*[0-9]+/) {
        split($2, a, /[ ,}]/)
        if (a[1] > 5) $2 = " 60" substr($2, length(a[1]) + 1)
    }
    print
    }
    ' $result_dir/treatment_alignments.filtered.json > ${result_dir}_mapq60/treatment_alignments.edited.json
    cp $result_dir/control_alignments.filtered.json ${result_dir}_mapq60
}

chipseq_graph() {
    local markname=$1
    local pipeline=$2
    local ref=$3
    local steps=$4
    local mapq_edit=$5
  
        # Steps of the pipeline
            # 1. Construct the graph
            # 2. Split the graph
            # 3. Align the reads
            # 4. Calculate parameters for peak calling
            # 5. Split json to chromosomes
            # 6. Call peaks
            # 7. Callpeaks_whole_genome_from_p_values + Find linear

    # Construct graphs
    if [[ " $steps " =~ 1 ]]; then
        variant_vcf=$wd/Graph_genome_data/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz
        backbone=$wd/Graph_genome_data/chm13v2.0.fa
        func="${pipeline}_graph"
        $func $variant_vcf $backbone $pipeline
    fi

    # Split the graph
    if [[ " $steps " =~ 2 ]]; then
        split_graph $pipeline $ref $wd/tools/gp.sif
    fi

    # Get input data
    input $markname
    
    if [ "$mapq_edit" == "mapq60" ]; then
        results_dir=$wd/results/${markname}/${pipeline}_${ref}_${mapq_edit}
        mkdir -p $results_dir
        edit_mapq $markname $pipeline $ref
        cp $wd/results/${markname}/${pipeline}_${ref}/param.txt $result_dir
        callpeaks "treatment_alignments.edited.json" "control_alignments.filtered.json" "$results_dir" "$ref" $pipeline "$steps" "narval" $wd/tools/gp.sif
    else     
        results_dir=$wd/results/${markname}/${pipeline}_${ref}
        # Alignment
        if [[ " $steps " =~ 3 ]]; then
            alignment $pipeline $ref $results_dir "$forward_trm" "$reverse_trm" "treatment"
            if [ -n "$forward_ctl" ]; then
                alignment $pipeline $ref $results_dir "$forward_ctl" "$reverse_ctl" "control"
            fi
        fi
        # Peak calling
        if [[ " $steps " =~ (^|[[:space:]])(4|5|6|7)($|[[:space:]]) ]]; then
            callpeaks "treatment_alignments.filtered.json" "control_alignments.filtered.json" "$results_dir" "$ref" $pipeline "$steps" "narval" $wd/tools/gp.sif
        fi
    fi
}