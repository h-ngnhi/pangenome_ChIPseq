#!/bin/bash 
#SBATCH --job-name=L1_map_146
#SBATCH --output=slurm-%j.out  # %j will be replaced with the job ID
#SBATCH --error=slurm-%j.err   # %j will be replaced with the job ID
#SBATCH --account=def-bourqueg
#SBATCH --time=8:00:00
#SBATCH --mem=249G
#SBATCH --cpus-per-task=64

# Enable xtrace to log the commands to stderr
set -x

start_time=$(date +%s)

# Set TMPDIR because the default /tmp is too small
export TMPDIR=/lustre07/scratch/hoangnhi/temp

export wd=/home/hoangnhi/projects/def-bourqueg/hoangnhi/ZNF146-507-Analysis-on-Pangenome
# graph_type=chm13
# data_dir=$wd/cgroza_data/H3K27AC_FLU
# graph_dir=$wd/Graph_genome_${graph_type}_trimmomatic
# results_dir=$data_dir/H3K27AC_${graph_type}_nocontrol
# vg_index=$graph_dir/giraffe_index

# mkdir -p $results_dir
# mkdir -p $results_dir/graphs
# mkdir -p $results_dir/graphs_pruned
# mkdir -p $results_dir/alignments
# mkdir -p $results_dir/vcf
# cd $results_dir
# ref=$wd/Graph_genome_data/chm13.draft_v1.1.fasta
# hprc_file=$wd/Graph_genome_data/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz
# vg autoindex --workflow giraffe -r $ref -v $hprc_file -p giraffe_index --threads 64
# vg autoindex --workflow giraffe -r $ref -p ${graph_type}_giraffe_index --threads 64
# module load samtools
# # Define the chromosomes
# chromosomes="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrM"
# IFS=',' read -r -a chromosome_array <<< "$chromosomes"
# # for chr in "${chromosome_array[@]}"; do
# #     # samtools faidx $ref $chr > graphs/${chr}.fa
# #     # vg construct -r graphs/${chr}.fa -m 32 -p > graphs/${chr}.vg
# #     # vg ids -j -m graphs/${chr}_mapping.txt graphs/${chr}.vg
# #     # vg prune -u -a -m graphs/${chr}_mapping.txt graphs/${chr}.vg -t 40 > graphs/${chr}.pruned.vg
# #     vg index -x graphs/${chr}.xg graphs/${chr}.vg -t 40
# #     vg index -k 16 -g graphs/${chr}.gcsa -p graphs/${chr}.xg graphs/${chr}.vg -t 40
# #     vg map -x graphs/${chr}.xg -g graphs/${chr}.gcsa -f $data_dir/Treatment.fastq.gz -t 40 -u 1 -m 1 > alignments/${chr}_alignments.gam
# #     vg surject -x graphs/${chr}.xg -b alignments/${chr}_alignments.gam > alignments/${chr}_alignments.bam
# # done
# module load tabix
# tabix -p vcf $hprc_file
# (seq 1 22; echo X; echo M) | parallel -j 24 "tabix -h $hprc_file chr{} > vcf/chr{}.vcf ; bgzip vcf/chr{}.vcf ; tabix vcf/chr{}.vcf.gz"
# (seq 1 22; echo X; echo M) | parallel -j 6 "vg construct -S -a -p -C -R chr{} -v vcf/chr{}.vcf.gz -r $ref -t 1 -m 32 > graphs/chr{}.vg"
# -----
#1. Construct the graph chm13 - no gbwt (haplotype-aware index) file because no variation
# (seq 1 22; echo X; echo M) | parallel -j 6 "vg construct -S -a -p -C -R chr{} -r $ref -t 1 -m 32 > graphs/chr{}.vg"
# vg ids -j -m graphs/mapping.txt graphs/chr*.vg
# (seq 1 22; echo X; echo M) | parallel -j 6 "vg prune -u -a -m graphs/mapping.txt graphs/chr{}.vg > graphs_pruned/chr{}.pruned.vg"
# vg index -g graphs/graph.gcsa -f graphs/mapping.txt graphs_pruned/*.vg -t 40
# vg map -x graphs/graph.xg -g graphs/graph.gcsa -f $data_dir/Treatment.fastq.gz -t 40 -u 1 -m 1 > alignments/treatment_alignments.gam
# vg surject -x graphs/graph.xg -b alignments/treatment_alignments.gam > alignments/treatment_alignments.bam


# (seq 1 22; echo X; echo M) | parallel -j 8 "touch -h vcf/chr{}.vcf.gz.tbi ; vg index -G graphs/chr{}.gbwt -v vcf/chr{}.vcf.gz graphs/chr{}.vg"
# vg gbwt -m -f -o graphs/vcfbub_graph.gbwt graphs/chr*.gbwt
# vg index -x graphs/vcfbub_graph.xg graphs/chr*.vg 
# (seq 1 22; echo X; echo M) | parallel -j 6 "vg prune -u -a -m graphs/mapping.txt -g graphs/chr{}.gbwt graphs/chr{}.vg > graphs_pruned/chr{}.pruned.vg"
# vg index -g graphs/vcfbub_graph.gcsa -f graphs/mapping.txt -Z 4096 graphs_pruned/*.vg
# vg map -x graphs/vcfbub_graph.xg -g graphs/vcfbub_graph.gcsa -f $data_dir/Treatment.fastq.gz -t 40 -u 1 -m 1 > alignments/treatment_alignments.gam
# vg surject -x graphs/vcfbub_graph.xg -b alignments/treatment_alignments.gam > alignments/treatment_alignments.bam

# vg autoindex -r $ref -p map_index --threads 64


# vg construct -r $ref -a > chm13_graph.vg 
# vg index -x chm13_graph.xg -g chm13_graph.gcsa chm13_graph.vg
# vg ids -j -m node_mapping.txt chm13_graph.vg
# vg prune -u -a -m node_mapping.txt chm13_graph.vg -t 40 > chm13_pruned.vg
# vg index -x chm13_graph.xg chm13_pruned.vg -t 40
# vg index -k 16 -g chm13_graph.gcsa -p chm13_graph.xg -X 3 -t 20
# vg index -j chm13_graph.lcp -g chm13_graph.gcsa chm13_graph.xg -t 20
# vg map -x chm13_graph.xg -g chm13_graph.gcsa -f $data_dir/Treatment.fastq.gz -t 40 -u 1 -m 1 > treatment_alignments.gam

# ------------------
# # #2. Align control and ChIP-seq file + filter alignments (based on mapq score)
# # # dir=$wd/GenPipes_sets/
# vg giraffe -Z $vg_index.gbz -m $vg_index.min -d $vg_index.dist -f $data_dir/treatment/H3K27AC.forward_treatment_1.fastq.gz -f $data_dir/treatment/H3K27AC.reverse_treatment_1.fastq.gz -t 64 > $results_dir/treatment_alignments.gam
# # # vg giraffe -Z $vg_index.gbz -m $vg_index.min -d $vg_index.dist -f $data_dir/Input.fastq.gz -t 64 > $results_dir/control_alignments.gam

# # # # cat 146Rep1_alignments.gam 146Rep2_alignments.gam > 146_alignments.gam
# vg filter $results_dir/treatment_alignments.gam -fu -q 1 -t 64 > $results_dir/treatment_alignments.filtered.gam
# # # vg filter $results_dir/control_alignments.gam -fu -q 30 -t 64 > $results_dir/control_alignments.filtered.gam

# # # # #3. Convert alignments file to json for Graph_peak_caller
# vg view -aj $results_dir/treatment_alignments.filtered.gam > $results_dir/treatment_alignments.filtered.json 
# # # vg view -aj $results_dir/control_alignments.filtered.gam > $results_dir/control_alignments.filtered.json

# # #4. Convert GAM to BAM and extract mapq scores for histogram
# vg surject -x $vg_index.gbz -b $results_dir/treatment_alignments.gam > $results_dir/treatment_alignments.bam

# module load samtools
# samtools view $results_dir/treatment_alignments.bam | awk '{print $5}' > $results_dir/mapq_scores.txt
# # awk '$1 ==0 {count++} END {print "MAPQ = 0:", count+0}' $results_dir/mapq_scores.txt
# # awk '$1 >0 && $1 < 30 {count++} END {print "0 < MAPQ < 30:", count+0}' $results_dir/mapq_scores.txt
# # awk '$1 >= 30 && $1 < 60 {count++} END {print "30 <= MAPQ < 60:", count+0}' $results_dir/mapq_scores.txt
# # awk '$1 == 60 {count++} END {print "MAPQ = 60:", count+0}' $results_dir/mapq_scores.txt


# # # # 5. Split the graph into chromosomes
# # # mkdir -p $results_dir/graphs
# # # vg chunk -C -x $vg_index.gbz -b $results_dir/graphs/

# # # ## Change the name of the graph files to the chromosome name
# # # for file in $results_dir/graphs/*.vg; do
# # #     chr_name=$(vg paths -L -x "$file" | head -n 1)
# # #     mv "$file" "graphs/${chr_name}.vg"
# # # done

# genome_size=3100000000
# # fragment_length=120 #calculated by 
# #     # SIF_IMAGE="$wd/macs2_latest.sif"
# #     # BIND_DIR="/lustre06/project/6002326/hoangnhi/ZNF146-507-Analysis-on-Pangenome:/mnt"
# #     # apptainer exec --contain --cleanenv --bind $BIND_DIR $SIF_IMAGE macs2 predictd -i /mnt/cgroza_data/H3K27AC/H3K27AC_CHM13graph/treatment_alignments.bam
# # read_length=102     #calculated by 
# # # read_length=$(cat "cgroza_data/H3K27AC/treatment/H3K27AC.forward_treatment_1.fastq.gz" | head -2 | tail -1 | wc -m)
# # unique_reads=13295105 #$(pcregrep --buffer-size 1000K -o1 '"sequence": "([ACGTNacgtn]{20,})"' cgroza_data/H3K27AC/H3K27AC_CHM13graph/treatment_alignments.filtered.json | sort | uniq | wc -l)
# # # echo $unique_reads
# # # unique_reads=112938149

# # # # Prepare the parameters file for using inside the container
# # # cat <<EOL > $results_dir/params.txt
# # # genome_size=$genome_size
# # # fragment_length=$fragment_length
# # # read_length=$read_length
# # # unique_reads=$unique_reads
# # # EOL
# module load bcftools
# bcftools view -i 'AC>1' $wd/Graph_genome_data/hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz -o $wd/Graph_genome_data/hprc-v1.1-mc-chm13.vcfbub.nosingleton.vcf.gz

vg_giraffe_graph() {            # To create graph and index for vg giraffe
    local variant_vcf=$1
    local ref=$2
    local graph_type=$3         # Name of the graph
    graph_dir=$wd/Pangenomes/vg_giraffe/$graph_type
    mkdir -p $graph_dir
    cd $graph_dir || exit
    if [ ! -f $variant_vcf ]; then
        vg autoindex --workflow giraffe -r $ref -p $graph_type --threads 64
    else
        vg autoindex --workflow giraffe -r $ref -v $variant_vcf -p $graph_type --threads 64
    fi
    cd $wd || exit
}
# variant_vcf=$wd/Graph_genome_data/L1_annotation/L1_annotation.vcfbub.vcf.gz
# ref=$wd/Graph_genome_data/chm13.draft_v1.1.fasta
# vg_giraffe_graph $variant_vcf $ref

vg_map_graph() {                 # To create graph and index for vg map
    local pipeline=$1           # vg_giraffe or vg_map
    local graph=$2              # chm13 or vcfbub

    ref=$wd/Graph_genome_data/chm13.draft_v1.1.fasta
    graph_dir=$wd/Pangenomes/$pipeline/$graph
    echo $graph_dir
    mkdir -p $graph_dir
    mkdir -p $graph_dir/graphs
    cd $graph_dir || exit
    if [ "$graph" == "chm13" ]; then
        (seq 1 22; echo X; echo M) | parallel -j 6 "vg construct -S -a -p -C -R chr{} -r $ref -t 1 -m 32 > graphs/chr{}.vg"
        vg ids -j -m mapping.txt graphs/chr*.vg
        vg index -x graph.xg graphs/chr*.vg -t 40
        vg index -g graph.gcsa -f mapping.txt graphs/chr*.vg -t 40
    else
        hprc_file=$wd/Graph_genome_data/L1_annotation/L1_annotation.vcfbub.vcf.gz
        module load tabix
        mkdir -p vcf
        (seq 1 22; echo X; echo M) | parallel -j 24 "tabix -h $hprc_file chr{} > vcf/chr{}.vcf ; bgzip vcf/chr{}.vcf ; tabix vcf/chr{}.vcf.gz"
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

# bgzip Graph_genome_data/L1_annotation/L1_annotation.vcfbub.vcf
# tabix -p vcf Graph_genome_data/L1_annotation/L1_annotation.vcfbub.vcf.gz

# vg_map_graph vg_map L1_vcfbub

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
alignment $pipeline $ref $results_dir "$forward" "$reverse" "treatment"
# 146_control
forward="$data_dir/146Input.trim.pair1.fastq.gz"
reverse="$data_dir/146Input.trim.pair2.fastq.gz"
alignment $pipeline $ref $results_dir "$forward" "$reverse" "control"

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