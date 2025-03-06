#!/bin/bash

#1. Split interleaved fastq file into forward and reverse fastq files
# If input fastq contain both forward and reverse
slit_interleavedfastq() {
    set=treatment
    treatment_set=("AF30_Flu_2-669476_ChIPmentation-H3K27ac_A00266_0289_1_S12_L001_I.fastq.gz")
    for num in $(seq 1 $num_treatment); do
        input_file="$dir/$set/${treatment_set[$num-1]}"
        if [[ -f $input_file ]]; then
            awk -v dir="$dir" -v set="$set" -v markname="$markname" -v num="$num" '
                NR%8==1 || NR%8==2 || NR%8==3 || NR%8==4 {print > dir "/" set "/" markname ".forward_" set "_" num ".fastq.gz"} 
                NR%8==5 || NR%8==6 || NR%8==7 || NR%8==0 {print > dir "/" set "/" markname ".reverse_" set "_" num ".fastq.gz"}
                ' "$input_file"
        else
            echo "File $input_file does not exist"
        fi  
    done

    set=control
    control_set=AF30_NI_2-669475_ChIPmentation-H3K27ac_A00266_0289_1_S11_L001_I.fastq.gz
    awk -v dir="$dir" -v set="$set" -v markname="$markname" '
        NR%8==1 || NR%8==2 || NR%8==3 || NR%8==4 {print > dir "/" set "/" markname ".forward_" set ".fastq"} 
        NR%8==5 || NR%8==6 || NR%8==7 || NR%8==0 {print > dir "/" set "/" markname ".reverse_" set ".fastq"}
        ' "$dir/$set/${control_set[0]}"
}

#2. readset file creator
create_readset() {
    local markname=$1
    local forward_trm=($2)
    local reverse_trm=($3)
    local forward_ctl=($4)
    local reverse_ctl=($5)
    local result_dir=$6

    num_treatment=${#forward_trm[@]}
    # Output TSV file name
    output_tsv="$result_dir/${markname}.readset.tsv"
    echo $output_tsv

    # Write the header of the TSV file
    echo -e "Sample\tReadset\tMarkName\tMarkType\tLibrary\tRunType\tRun\tLane\tAdapter1\tAdapter2\tQualityOffset\tBED\tFASTQ1\tFASTQ2\tBAM" > $output_tsv

    # Write the treatment line(s) to the TSV file
    for num in $(seq 1 $num_treatment); do
        sample="${markname}_treatment_$num"
        readset="treatment_$num"
        fastq1="${forward_trm[$((num-1))]}.gz"
        fastq2="${reverse_trm[$((num-1))]}.gz"

        echo -e "$sample\t$readset\t$markname\tN\t\tPAIRED_END\trun10\t1\t\t\t33\t\t$fastq1\t$fastq2\t" >> $output_tsv
    done

    # Write the control line to the TSV file
    sample="${markname}_control"
    readset="control"
    fastq1="${forward_ctl[0]}.gz"
    fastq2="${reverse_ctl[0]}.gz"

    echo -e "$sample\t$readset\tInput\tI\t\tPAIRED_END\trun10\t1\t\t\t33\t\t$fastq1\t$fastq2\t" >> $output_tsv

    # Notify user
    echo "TSV file '$output_tsv' created successfully."
}

#3. design file creator
create_design() {
    local markname=$1
    local result_dir=$2

    # Output file name
    output_txt="$result_dir/${markname}.design.txt"

    # Write the header to the design file
    echo -e "Sample\tMarkName\t${markname}_vs_${markname}Input" > $output_txt

    # Write the treatment lines to the design file
    for num in $(seq 1 $num_treatment); do
        sample="$markname"
        echo -e "$sample\t$markname\t0" >> $output_txt
    done

    # Write the control line to the design file
    sample="$markname"
    marktype="Input"
    echo -e "$sample\t$marktype\t0" >> $output_txt

    # Notify user
    echo "Design file '$output_txt' created successfully."
}

#4. config file creator
# create_genome() {
#     local ref=$1            # Directory of FASTA reference

#     module load samtools
#     module load bwa

#     # Create genome dictionary file
#     samtools dict $ref -o ${ref%.*}.dict

#     # Create BWA index
#     bwa index $ref

#     # Create chromosome size file
#     samtools faidx $ref -o ${ref%.*}.dict
# }
create_config() {
    local ref=$1                    # Directory for linear reference FASTA
    local result_dir=$2
    output_ini_file="$result_dir/t2t.ini"

    # Create and write content to the ini file
    cat <<EOL > $output_ini_file
[DEFAULT]
assembly=T2T-CHM13v2.0
genome_fasta=$ref
genome_dictionary=${ref%.*}.dict
genome_bwa_index=${ref%.*}.bwaindex.fa
chromosome_size=$ref.fai

[trimmomatic]
min_length = 98

[sambamba_view_filter]
min_mapq = 30

[macs2_callpeak]
extsize = 150
other_options = --keep-dup all
EOL

    # Notify user
    echo "INI file '$output_ini_file' created successfully."
}

#5. Run chipseq
chipseq_gp() {
    local markname=$1
    local result_dir=$2

    cd $results_dir
    #$MUGQIC_PIPELINES_HOME only available to ls command after loading mugqic/genpipes/4.4.5
    genpipes chipseq -c $wd/tools/genpipes/genpipes/pipelines/chipseq/chipseq.base.ini \
                $wd/genome_data/GenPipes_t2t/narval.ini $result_dir/t2t.ini \
                -r $result_dir/${markname}.readset.tsv -d $result_dir/${markname}.design.txt \
                -o $result_dir > chipseqScript.sh
    bash chipseqScript.sh
    cd $wd
}
# Check if the trimming adapter is needed - normally no for long dna fragments
# head -n 10000 "$dir/$set/forward_$set.fastq" | grep -E 'AGATCGGAAGAGC'
# for treatment
