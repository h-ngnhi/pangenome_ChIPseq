include { create_readset_tsv; create_design_txt; create_config_ini; chipseq_linear }  from './modules/linear_Genpipes.nf'


process read_input_file {
    input:
    path input_file

    output:
    tuple val(markname), val(library), val(samples), val(fastq1s), val(fastq2s)

    script:
    def lines = input_file.toFile().readLines()
    def markname = lines[1].split("\t")[1]
    def library = lines[1].split("\t")[2]
    def samples = lines.drop(1).collect { it.split("\t")[0] }
    def fastq1s = lines.drop(1).collect { it.split("\t")[3] }
    def fastq2s = lines.drop(1).collect { it.split("\t")[4] }
}

// Display help message if requested
if (params.help || !params.input || !params.ref) {
    log.info """
    Usage: nextflow run hoangnhi-nguyen/pangenome_ChIPseq [options]

    Options:
    --input             Path to the input tsv file (required)
                        The input file should have the following columns:
                        Sample, MarkName, Library, FASTQ1, FASTQ2
                        Paths to the FASTQ files should be relative to the input file 
    --ref               Reference graph type ('graph' or 'linear', default: 'graph')

    --cpus              Number of CPUs to allocate (default: '64')
    --memory            Memory allocation (default: '249 GB')
    --time              Time allocation (default: '24h')

    --help              Show this help message and exit

    Example:
    nextflow run h-ngnhi/pangenome_ChIPseq --input input.txt --ref linear --cpus 64 --memory '249 GB' --time '24h'
    """
    exit 0
}

workflow {
    Channel.fromPath(params.input, checkIfExists:true).set{input_file_ch}
    def read_input_file_out = read_input_file(input_file_ch)
    log.info "Account parameter: ${params.account}"
    if (params.ref=="linear") {
        log.info "Running linear ChIP-Seq pipeline"
        def markname_ch = read_input_file_out.map { it[0] }
        def readset_ch = create_readset_tsv(read_input_file_out)
        def design_ch = create_design_txt(read_input_file_out)
        def ini_ch = create_config_ini()
        readset_ch.combine(markname_ch).combine(design_ch).combine(ini_ch).set{chipseq_input}
        chipseq_linear(chipseq_input)
    }
}
