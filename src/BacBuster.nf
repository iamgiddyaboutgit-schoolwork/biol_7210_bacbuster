#!/usr/bin/env nextflow
// In the future, it looks like the default will be DSL2,
// so we may as well get used to it.
nextflow.enable.dsl=2

params.seq_reads = "/home/team3/raw_data/Raw_FQs/*"

raw_reads = Channel.fromPath(params.seq_reads, checkIfExists: true)

process pre_assembly_reads_qc {
    conda "Team3-WebServer_falco_env.yml"
    input:
    // https://www.nextflow.io/docs/edge/process.html#input-type-file
    path raw_fastq

    output:
    file "${raw_fastq.baseName}.txt"

    // https://www.nextflow.io/docs/edge/process.html#shell
    shell:
    '''
    falco --outdir ./ --threads 1 -subsample 1000 -skip-data -skip-report !{raw_fastq}
    '''
}

workflow {
    pre_assembly_reads_qc(raw_reads)
}