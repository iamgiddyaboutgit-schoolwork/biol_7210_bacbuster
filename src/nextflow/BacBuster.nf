#!/usr/bin/env nextflow
// In the future, it looks like the default will be DSL2,
// so we may as well get used to it.
nextflow.enable.dsl=2

process trim {
    conda "Team3-WebServer_env.yml"
    // https://www.nextflow.io/docs/edge/process.html#input-type-file
    input:
        // https://www.nextflow.io/docs/latest/channel.html#fromfilepairs
        tuple val(sample_id), path("fq")
    output:
        tuple path("${sample_id}_trimmed_1.fq.gz"), path("${sample_id}_trimmed_2.fq.gz")
    shell:
    """
    #!/usr/bin/env bash
 
    # If no quality encoding is specified,
    # it will be determined automatically (since version 0.32). The prior default was -phred64.
    # https://github.com/cemalley/Ryu_cerebral_organoids/blob/96543e0e835b1c8a13bacb45f2b34fa47f1cf3dc/Bulk_RNASeq/Swarm_single_trim.sh
    # Bundled Illumina adapter FASTA files can be found at the path where trimmomatic is installed, such as:
    # mambaforge/envs/Team3-WebServer_env/share/trimmomatic-0.39-2/adapters
    # You can determine which adapters were used by running and comparing the output.
    # If the correct adapter FASTA is chosen, then more trimming should be performed.
    # More information on Illumina truseq_dna_sample_prep_kits dated 2014:
    # https://www.illumina.com/documents/products/datasheets/datasheet_truseq_dna_sample_prep_kits.pdf
    # We probably do not have TruSight adapter sequences because they seem to be only marketed for 
    # human samples: https://www.illumina.com/products/trusight-panels.html
    # We most likely have TruSeq adapters. More info. on Illumina adapter sequences
    # can be found in their PDF here: https://support.illumina.com/downloads/illumina-adapter-sequences-document-1000000002694.html
    # Make the MAXINFO targetLength a little higher than the recommended 40 to 
    # more heavily penalize reads that are shorter than targetLength.
    # Note that Illumina reads should all have the same length,
    # so this parameter is likely irrelevant except when trimming is done that shortens the read lengths.
    # https://www.ecseq.com/support/ngs/why-do-the-reads-all-have-the-same-length-when-sequencing-differently-sized-fragments#:~:text=So%20why%20you%20get%20reads,1).
    # strictness is on a scale from 0 to 1, with values closer to 1 causing more aggressive trimming.

    trimmomatic PE \
        -threads 4 \
        !{fq[0]} \
        !{fq[1]} \
        !{sample_id}_trimmed_1.fq.gz \
        !{sample_id}_trimmed_1_unpaired.fq.gz \
        !{sample_id}_trimmed_2.fq.gz \
        !{sample_id}_trimmed_2_unpaired.fq.gz \
        ILLUMINACLIP:/home/jpatterson87/bin/mambaforge/envs/Team3-WebServer_env/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa \
        MAXINFO:50:0.50
    """
}

workflow {
    params.seq_reads = "/home/jpatterson87/big_project/Team3-WebServer/testing_data/sequencing_reads/*_{1,2}.fq.gz"
    // https://www.nextflow.io/docs/latest/channel.html#fromfilepairs
    // raw_reads = Channel.fromFilePairs(params.seq_reads, maxDepth=1, checkIfExists: true)
    raw_reads = Channel.fromFilePairs(params.seq_reads, maxDepth:1, checkIfExists:true)
    trim(raw_reads) 
}