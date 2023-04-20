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
        tuple path("${sample_id}.trimmed_1.fq.gz"), path("${sample_id}.trimmed_2.fq.gz")
    shell:
    """
    #!/usr/bin/env bash
 
    # More info. on Illumina adapter sequences
    # can be found in their PDF here: https://support.illumina.com/downloads/illumina-adapter-sequences-document-1000000002694.html

    # sample_id will be like CGT3027.downsampled according to the glob pattern.

    # Assume reads are paired end.
    # Assume quality encoding is phred33.  TODO: check this.

    # fastp: https://doi.org/10.1093/bioinformatics/bty560

    # Low complexity regions in the proteins of 
    # prokaryotes perform important functional roles and are highly conserved:
    # https://doi.org/10.1093/nar/gkz730

    # Filtering by average quality happens after all trimming according to:
    # https://github.com/OpenGene/fastp/issues/159

    fastp \
        --in1 !{fq[0]} \
        --out1 !{sample_id}.trimmed_1.fq.gz \
        --in2 !{fq[1]} \
        --out2 !{sample_id}.trimmed_2.fq.gz \
        --detect_adapter_for_pe \
        --disable_trim_poly_g \
        --cut_right \
        --cut_window_size 4 \
        --cut_mean_quality 25 \
        --length_required 15 \
        --qualified_quality_phred 17 \
        --unqualified_percent_limit 30 \
        --n_base_limit 20 \
        --thread 4
    """
}

workflow {
    params.seq_reads = "/home/jpatterson87/big_project/Team3-WebServer/testing_data/sequencing_reads/*_{1,2}.fq.gz"
    // https://www.nextflow.io/docs/latest/channel.html#fromfilepairs
    // raw_reads = Channel.fromFilePairs(params.seq_reads, maxDepth=1, checkIfExists: true)
    raw_reads = Channel.fromFilePairs(params.seq_reads, maxDepth:1, checkIfExists:true)
    trim(raw_reads) 
}