include * from './linear_Genpipes.nf'

process read_input_file {
    input:
    path input_file

    output:
    tuple val(markname), val(library), val(samples), val(fastq1s), val(fastq2s)

    script:
    """
    # Extract the values from the input file
    markname=\$(awk 'NR==2 {print \$2}' ${input_file})   # Get the MarkName from the second column (assumed consistent)
    library=\$(awk 'NR==2 {print \$3}' ${input_file})    # Get the Library from the third column (assumed consistent)

    # Extract the list of Samples, FASTQ1, and FASTQ2 files
    samples=(\$(awk 'NR>1 {print \$1}' ${input_file} | tr '\\n' ' '))
    fastq1s=(\$(awk 'NR>1 {print \$4}' ${input_file} | tr '\\n' ' '))
    fastq2s=(\$(awk 'NR>1 {print \$5}' ${input_file} | tr '\\n' ' '))
    """
}
// Ensure that users provide required parameters
if (!params.markname || params.treatment_fastq_fw.size() == 0 || params.treatment_fastq_rv.size() == 0 || params.control_fastq_fw.size() == 0 || params.control_fastq_rv.size() == 0) {
    exit 1, "ERROR: Required parameters are missing. Please provide 'markname', 'treatment_fastq_fw', 'treatment_fastq_rv', 'control_fastq_fw', and 'control_fastq_rv'."
}

workflow {
    Channel.fromPath(params.input, checkIfExists:true).set{input_file_ch}
    read_input_file(input_file_ch)
    if (params.ref=="linear") {
        log.info "Running linear ChIP-Seq pipeline"
        create_readset_tsv(read_input_file.out)
        create_design_txt(read_input_file.out)
        create_config_ini()
        chipseq_linear(create_readset_tsv.out, create_design_txt.out, create_config_ini.out)
        log.info "Linear ChIP-Seq pipeline completed"
    }
}