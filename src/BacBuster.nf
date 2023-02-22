#!/usr/bin/env nextflow
// In the future, it looks like the default will be DSL2,
// so we may as well get used to it.
nextflow.enable.dsl=2

params.seq_reads = "/home/team3/raw_data/Raw_FQs/*"

raw_reads = Channel.fromPath(params.seq_reads)

process pre_assembly_reads_qc {
    mamba "Team3-WebServer_falco_env.yml"
    input:
    file raw_fastq from raw_reads

    output:
    file "${raw_fastq.baseName}.txt" into raw_reads_qc

    shell:
    """
    falco --outdir ./ --threads 1 -subsample 1000 ${raw_fastq}
    """
}

workflow {
    = pre_assembly_reads_qc(raw_reads)
}