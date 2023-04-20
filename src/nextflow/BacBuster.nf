#!/usr/bin/env nextflow
// In the future, it looks like the default will be DSL2,
// so we may as well get used to it.
nextflow.enable.dsl=2

params.seq_reads = "/home/jpatterson87/big_project/Team3-WebServer/testing_data/sequencing_reads/*_{1,2}.fq.gz"

process trim {
    conda "Team3-WebServer_env.yml"
    // https://www.nextflow.io/docs/edge/process.html#input-type-file
    input:
        tuple val(sample_id), path(fq)
    output:
        stdout
    shell:
    """
    #!/usr/bin/env bash

    # Assume that 
    # trimmomatic PE \
    #     -threads 4 \
    
    """
}

workflow {
    raw_reads = Channel.fromFilePairs(params.seq_reads, checkIfExists: true)
    trim(raw_reads) | view
}