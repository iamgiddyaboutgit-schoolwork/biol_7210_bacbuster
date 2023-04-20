#!/usr/bin/env nextflow
// In the future, it looks like the default will be DSL2,
// so we may as well get used to it.
nextflow.enable.dsl=2

params.seq_reads = "/home/jpatterson87/big_project/Team3-WebServer/testing_data/sequencing_reads/*.fq.gz"

raw_reads = Channel.fromPath(params.seq_reads, checkIfExists: true)



workflow {
    
}