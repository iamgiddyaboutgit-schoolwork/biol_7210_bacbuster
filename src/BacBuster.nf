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
    // The .baseName part here is used to get the file name 
    // without its extension e.g. /some/path/file.tar.gz -> file.tar
    file "${raw_fastq.baseName}.falco_summary.txt"

    // https://www.nextflow.io/docs/edge/process.html#shell
    shell:
    '''
    # -summary-filename: The isolate name is extracted for naming purposes 
    falco \
        --outdir ./ \
        --threads 1 \
        -subsample 1000 \
        -skip-data \
        -skip-report \
        -summary-filename $(basename !{raw_fastq} | cut --delimiter=. --fields=1).falco_summary.txt \
        !{raw_fastq}
    '''
}

workflow {
    pre_assembly_reads_qc(raw_reads)
}