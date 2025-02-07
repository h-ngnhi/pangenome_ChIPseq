#!/bin/bash 

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
        (seq 1 22; echo X; echo M) | parallel -j 24 "tabix -h $variant_vcf chr{} > vcf/chr{}.vcf ; bgzip vcf/chr{}.vcf ; tabix vcf/chr{}.vcf.gz"
        (seq 1 22; echo X; echo M) | parallel -j 6 "vg construct -S -a -p -C -R chr{} -v vcf/chr{}.vcf.gz -r $ref -t 1 -m 32 > graphs/chr{}.vg"
        vg ids -j -m mapping.txt graphs/chr*.vg
        # Adding variants --> need to create gbwt (haplotype-aware index) file
        (seq 1 22; echo X; echo M) | parallel -j 8 "touch -h vcf/chr{}.vcf.gz.tbi ; vg index -G graphs/chr{}.gbwt -v vcf/chr{}.vcf.gz graphs/chr{}.vg"
        vg gbwt -m -f -o graph.gbwt graphs/chr*.gbwt
        vg index -x graph.xg graphs/chr*.vg 
        # Prune the graph for gcas index
        mkdir -p graphs_pruned
        (seq 1 22; echo X; echo M) | parallel -j 6 "vg prune -u -a -m mapping.txt -g graphs/chr{}.gbwt graphs/chr{}.vg > graphs_pruned/chr{}.pruned.vg"
        vg index -g graph.gcsa -f mapping.txt -Z 4096 graphs_pruned/*.vg -t 40
    fi
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
    vg filter ${input}_alignments.gam -fu -q 30 -t 64 > ${input}_alignments.filtered.gam
    vg view -aj ${input}_alignments.filtered.gam > ${input}_alignments.filtered.json
    if [ "$input" == "treatment" ]; then
        module load samtools
        samtools view ${input}_alignment.bam | awk '{print $5}' > mapq_scores.txt
        awk '$1 ==0 {count++} END {print "MAPQ = 0:", count+0}' mapq_scores.txt
        awk '$1 >0 && $1 < 30 {count++} END {print "0 < MAPQ < 30:", count+0}' mapq_scores.txt
        awk '$1 >= 30 && $1 < 60 {count++} END {print "30 <= MAPQ < 60:", count+0}' mapq_scores.txt
        awk '$1 == 60 {count++} END {print "MAPQ = 60:", count+0}' mapq_scores.txt
    fi
    cd $wd || exit
}

# samtools view 146_vgmap_results/L1_vcfbub/treatment_alignments.bam | awk '{print $5}' > mapq_scores.txt
# awk '$1 ==0 {count++} END {print "MAPQ = 0:", count+0}' mapq_scores.txt
# awk '$1 >0 && $1 < 30 {count++} END {print "0 < MAPQ < 30:", count+0}' mapq_scores.txt
# awk '$1 >= 30 && $1 < 60 {count++} END {print "30 <= MAPQ < 60:", count+0}' mapq_scores.txt
# awk '$1 == 60 {count++} END {print "MAPQ = 60:", count+0}' mapq_scores.txt
pipeline=vg_giraffe
ref=L1_vcfbub
data_dir=$wd/Graph_genome_data
results_dir=$wd/146_giraffe_results/$ref
# data_dir=cgroza_data/H3K27AC_FLU
# 146_treatment
forward="$data_dir/146Rep1.trim.pair1.fastq.gz $data_dir/146Rep2.trim.pair1.fastq.gz"
reverse="$data_dir/146Rep1.trim.pair2.fastq.gz $data_dir/146Rep2.trim.pair2.fastq.gz"
# forward="$wd/$data_dir/treatment/H3K27AC.forward_treatment_1.fastq.gz"
# reverse="$wd/$data_dir/treatment/H3K27AC.reverse_treatment_1.fastq.gz"
# alignment $pipeline $ref $results_dir "$forward" "$reverse" "treatment"
# 146_control
forward="$data_dir/146Input.trim.pair1.fastq.gz"
reverse="$data_dir/146Input.trim.pair2.fastq.gz"
# alignment $pipeline $ref $results_dir "$forward" "$reverse" "control"

split_graph() {                         # To split graphs by chromosomes for GP
    local pipeline=$1                   # vg_giraffe or vg_map
    local graph=$2                      # chm13 or vcfbub or hprc-v1.1-mc-chm13 or L1_vcfbub

    graph_dir=$wd/Pangenomes/$pipeline/$graph
    if [ "$pipeline" == "vg_giraffe" ]; then
        #Split the graph into chromosomes
        mkdir -p $graph_dir/graphs
        vg chunk -C -x $graph_dir/$graph.gbz -b $graph_dir/graphs/

        # Change the name of the graph files to the chromosome name
        for file in $graph_dir/graphs/*.vg; do
            chr_name=$(vg paths -L -x $file | head -n 1)
            mv $file $graph_dir/graphs/${chr_name}.vg
        done
    fi

    graph_dir=$graph_dir/graphs
    Convert graphs to json
    (seq 1 22; echo X; echo M) | parallel -j 2 "vg view -j $graph_dir/chr{}.vg > $graph_dir/chr{}.json ; vg stats -r $graph_dir/chr{}.vg | cut -f 2 > $graph_dir/node_range_chr{}.txt"

    # Create ob_graph files
    SIF_IMAGE="$wd/gp.sif"
    BIND_DIR="/lustre06/project/6002326/hoangnhi/ZNF146-507-Analysis-on-Pangenome/${graph_dir#$wd/}:/mnt_data"

    module load apptainer
    apptainer exec --contain --cleanenv --bind $BIND_DIR $SIF_IMAGE /bin/bash -c "
        for file in /mnt_data/*.json; do
            graph_peak_caller create_ob_graph \$file
        done
    "
}
# split_graph $pipeline $ref

callpeaks() {
    local treatment=$1              # treatment json file
    local control=$2                # control json file (could be empty if no control)
    local results_dir=$3            # result directory
    local graph=$4                  # chm13 or vcfbub
    local pipeline=$5               # vg_giraffe or vg_map
    local fragment_length=$6        # fragment length
    local read_length=$7            # read length
    local unique_reads=$8           # unique reads

    # Prepare the parameters file for using inside the container (no indent for this command)
cat <<EOL > $wd/params.txt
fragment_length=$fragment_length
read_length=$read_length
unique_reads=$unique_reads
graph=$pipeline/$graph
results_dir=$results_dir
treatment=$treatment
control=$control
EOL

    SIF_IMAGE="$wd/gp.sif"
    BIND_DIR="/lustre06/project/6002326/hoangnhi/ZNF146-507-Analysis-on-Pangenome:/mnt_data"

    module load apptainer
    apptainer exec --contain --cleanenv --bind $BIND_DIR $SIF_IMAGE /bin/bash -c "
        while IFS= read -r line; do 
            export \$line
        done < /mnt_data/params.txt
        results_dir=/mnt_data/\$results_dir
        graph_dir=/mnt_data/Pangenomes/\$graph/graphs

        # Split the treatment and control json files into chromosomes
        # at this step, it is required to pass in all the chromosomes at ones to the command
        chromosomes=\"chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrM\"
        graph_peak_caller split_vg_json_reads_into_chromosomes \$chromosomes \$results_dir/\$treatment \$graph_dir/
        if [ -n \"\$control\" ]; then
            graph_peak_caller split_vg_json_reads_into_chromosomes \$chromosomes \$results_dir/\$control \$graph_dir/
        fi
        mkdir -p \$results_dir/reads_by_chrom
        mv \$results_dir/*chr*.json \$results_dir/reads_by_chrom/


        # Call peaks
        mkdir -p \$results_dir/callpeaks/
        chromosomes=\$(echo \$chromosomes | tr ',' ' ')
        for chromosome in \$chromosomes; do
            echo \$chromosome
            treatment_chr=\$(basename \${treatment%.*})
            echo \$treatment_chr
            if [ -n \"\$contr\" ]; then
                ctl_chr=\$(basename \${control%.*})
                echo if cond is reached \$ctl_chr
                graph_peak_caller callpeaks -g \$graph_dir/\$chromosome.nobg -s \$results_dir/reads_by_chrom/\${treatment_chr}_\$chromosome.json -c \$results_dir/reads_by_chrom/\${ctl_chr}_\$chromosome.json -f \$fragment_length -r \$read_length -u \$unique_reads -p True -G 3100000000 -n \$results_dir/callpeaks/\${chromosome}_
            else
                graph_peak_caller callpeaks -g \$graph_dir/\$chromosome.nobg -s \$results_dir/reads_by_chrom/\${treatment_chr}_\$chromosome.json -f \$fragment_length -r \$read_length -u \$unique_reads -p True -G 3100000000 -n \$results_dir/callpeaks/\${chromosome}_
            fi
        done

        cd \$results_dir/callpeaks/
        for chromosome in \$chromosomes; do
            echo \$chromosome
            graph_peak_caller callpeaks_whole_genome_from_p_values -d \$graph_dir/ -n \"\" -f \$fragment_length -r \$read_length \$chromosome
        done
        chromosomes=\"chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrM\"
        graph_peak_caller concatenate_sequence_files \$chromosomes all_peaks.fasta
        chromosomes=\$(echo \$chromosomes | tr ',' ' ') 
        for chromosome in \$chromosomes; do
            echo \$chromosome
            graph_peak_caller find_linear_path -g \$graph_dir/\$chromosome.nobg \$graph_dir/\$chromosome.json \$chromosome \$results_dir/reads_by_chrom/\${chromosome}_linear_path.interval
            graph_peak_caller peaks_to_linear \${chromosome}_max_paths.intervalcollection \$results_dir/reads_by_chrom/\${chromosome}_linear_path.interval \$chromosome \${chromosome}_linear_peaks.bed
        done
        mv \$graph_dir/.interval \$results_dir/reads_by_chrom/
        mv \$graph_dir/*.interval \$results_dir/reads_by_chrom/
        mv \$graph_dir/*.npz \$results_dir/reads_by_chrom/
    "
}