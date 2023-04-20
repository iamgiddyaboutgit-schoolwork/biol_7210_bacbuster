#!/usr/bin/env nextflow
// In the future, it looks like the default will be DSL2,
// so we may as well get used to it.
nextflow.enable.dsl=2

params.seq_reads = "/home/jpatterson87/big_project/Team3-WebServer/testing_data/sequencing_reads/*.fq.gz"

raw_reads = Channel.fromPath(params.seq_reads, checkIfExists: true)

process pre_assembly_reads_qc {
    conda "Team3-WebServer_env.yml"
    input:
    // https://www.nextflow.io/docs/edge/process.html#input-type-file
    path raw_fastq

    output:
    // The .baseName part here is used to get the file name 
    // without its extension e.g. /some/path/file.tar.gz -> file.tar
    file "summary.txt"

    // https://www.nextflow.io/docs/edge/process.html#shell
    shell:
    '''
    falco \
        --outdir . \
        --threads 1 \
        -subsample 1000 \
        -skip-data \
        -skip-report \
        !{raw_fastq}
    '''
}

workflow {
    pre_assembly_reads_qc(raw_reads)
}