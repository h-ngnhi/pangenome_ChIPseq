process create_readset_tsv {
    input:
    tuple val(markname), val(library), val(samples), val(fastq1s), val(fastq2s)

    output:
    path "${markname}.readset.tsv"

    script:
    """
    output_tsv="${markname}.readset.tsv"
    echo -e "Sample\tReadset\tMarkName\tMarkType\tLibrary\tRunType\tRun\tLane\tAdapter1\tAdapter2\tQualityOffset\tBED\tFASTQ1\tFASTQ2\tBAM" > \$output_tsv
    for i in \$(seq 0 \$((${#samples[@]}-1))); do
        if [[ \$sample == *"control"* ]]; then
            marktype="Input"
        else
            marktype=\${markname}
        fi
        echo -e "\${samples[\$i]}\t\${samples[\$i]}\t\${markname}\t\${marktype}\t\${library}\tPAIRED_END\trun10\t1\t\t\t33\t\t\${fastq1s[\$i]}\t\${fastq2s[\$i]}\t" >> \$output_tsv
    done
    """
}

process create_design_txt {
    input:
    tuple val(markname), val(samples)

    output:
    path "${markname}.design.txt"

    script:
    """
    output_txt="${markname}.design.txt"
    echo -e "Sample\tMarkName\t${markname}_vs_${markname}Input" > \$output_txt
    for i in \$(seq 0 \$((${#samples[@]}-2))); do
        sample=\${samples[\$i]}
        echo -e "\${samples[\$i]}\t\$markname\t0" >> \$output_txt
    done
    echo -e "\${samples[-1]}\tInput\t0" >> \$output_txt
    """
}

process create_config_ini {
    output:
    path "T2T.ini"

    script:
    """
    output_ini_file="T2T.ini"

    cat <<EOL > \$output_ini_file
[DEFAULT]
assembly=T2T-CHM13v2.0
assembly_dir=\$MUGQIC_INSTALL_HOME_DEV/genomes/species/%(scientific_name)s.%(assembly)s
genome_fasta=%(assembly_dir)s/genome/Homo_sapiens.T2T-CHM13v2.0.maskedY.rCRS.EBV.fa
genome_dictionary=%(assembly_dir)s/genome/Homo_sapiens.T2T-CHM13v2.0.maskedY.rCRS.EBV.dict
genome_bwa_index=%(assembly_dir)s/genome/bwa-mem_index/Homo_sapiens.T2T-CHM13v2.0.maskedY.rCRS.EBV.fa
chromosome_size=%(assembly_dir)s/genome/Homo_sapiens.T2T-CHM13v2.0.maskedY.rCRS.EBV.fa.fai

[trimmomatic]
min_length = 98

[sambamba_view_filter]
min_mapq = 30

[macs2_callpeak]
extsize = 150
other_options = --keep-dup all
EOL
    """
}

//Will modify when adding chipseq for graph
process chipseq_linear() {
    input:
    tuple path(readset_tsv), path(design_txt), path(ini_file)

    output:
    path "${params.markname}_CHM13linear"

    script:
    """
    module load mugqic/genpipes/4.4.5
    module load mugqic/python/3.10.4
    directory="${params.markname}_CHM13linear"
    mkdir -p \$directory
    #\$MUGQIC_PIPELINES_HOME only available to ls command after loading mugqic/genpipes/4.4.5
    chipseq.py -c \$MUGQIC_PIPELINES_HOME/pipelines/chipseq/chipseq.base.ini \
                \$MUGQIC_PIPELINES_HOME/pipelines/common_ini/narval.ini $ini_file \
                -r $readset_tsv -d $design_txt -o $directory > "chipseqScript_${params.markname}.txt" 
    bash "chipseqScript_${params.markname}.txt" 
    """
}