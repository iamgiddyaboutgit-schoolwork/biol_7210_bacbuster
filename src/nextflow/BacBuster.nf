#!/usr/bin/env nextflow
// In the future, it looks like the default will be DSL2,
// so we may as well get used to it.
nextflow.enable.dsl=2

process trim {
    conda "/home/andy/compGen/Team3-WebServer/src/nextflow/Team3-WebServer_env.yml"
    // https://www.nextflow.io/docs/edge/process.html#input-type-file
    input:
        // https://www.nextflow.io/docs/latest/channel.html#fromfilepairs
        tuple val(sample_id), path("fq")
    output:
        tuple val(sample_id), path("${sample_id}.trimmed_1.fq.gz"), path("${sample_id}.trimmed_2.fq.gz")
    shell:
    '''
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

    # Order of filtering:
    # https://github.com/OpenGene/fastp/

    fastp \
        --in1 !{fq[0]} \
        --out1 !{sample_id}.trimmed_1.fq.gz \
        --in2 !{fq[1]} \
        --out2 !{sample_id}.trimmed_2.fq.gz \
        --detect_adapter_for_pe \
        --disable_trim_poly_g \
        --unqualified_percent_limit 30 \
        --n_base_limit 20 \
        --html "!{sample_id}.fastp.html" \
        --report_title "!{sample_id} fastp Report" \
        --thread 4
    '''
}

process assemble {
    conda "/home/andy/compGen/Team3-WebServer/src/nextflow/skesa.yml"
    input:
        tuple val(sample_id), path("${sample_id}.trimmed_1.fq.gz"), path("${sample_id}.trimmed_2.fq.gz")
    output:
        tuple val(sample_id), path("${sample_id}.trimmed.contigs.fq.gz")
    shell:
    '''
    #!/usr/bin/env bash

    # SKESA is designed for ILLUMINA reads only.
    #
    # min_count
    # "All k-mers with count below the minimum count are ignored in the assembly."
    # A higher value for min_count appears to make the bloom filter use
    # more memory.  
    #
    # max_kmer
    # The default value for max_kmer appears to be 0.
    #
    # max_kmer_count
    # If max_kmer_count is not specified, then it appears to be estimated
    # from a starting value of 10.
    #
    # vector_percent
    # vector_percent appears to refer to the percentage of reads 
    # containing a given 19-mer for the 19-mer to be considered a vector.
    #
    # Skesa actually attempts to do adapter trimming before assembly.
    #
    # "No explicit error correction of reads is done by SKESA
    # as the heuristics of SKESA can handle the errors in a typical illumina read set."

    skesa \
        --cores 4 \
        --memory 8 \
        --reads !{sample_id}.trimmed_1.fq.gz,!{sample_id}.trimmed_2.fq.gz \
        --contigs_out !{sample_id}.trimmed.contigs.fq.gz \
        --steps 15 
    '''
}

process predict_genes {
    conda "/home/andy/compGen/Team3-WebServer/src/nextflow/Team3-WebServer_env.yml"
    // https://www.nextflow.io/docs/edge/process.html#input-type-file
    input:
        // https://www.nextflow.io/docs/latest/channel.html#fromfilepairs
        tuple val(sample_id), path("draft_assembly")
    output:
        tuple val(sample_id), path("${sample_id}.faa")
        tuple val(sample_id), path("${sample_id}.fna")
        tuple val(sample_id), path("${sample_id}.prodigal.out")
    shell:
    """
    #!/usr/bin/env bash
    
    prodigal \
        -a !{sample_id}.faa \
        -d !{sample_id}.fna \
        -f gff \
        -i !{draft_assembly} \
        -o !{sample_id}.prodigal.out

    """
}

workflow {
    //params.seq_reads = "/home/jpatterson87/big_project/Team3-WebServer/testing_data/sequencing_reads/*_{1,2}.fq.gz"
    //path = params.seq_reads
    params.seq_reads = ""
    // https://www.nextflow.io/docs/latest/channel.html#fromfilepairs
    // raw_reads = Channel.fromFilePairs(, maxDepth=1, checkIfExists: true)
    raw_reads = Channel.fromFilePairs(params.seq_reads, maxDepth:1, checkIfExists:true)
    trim(raw_reads) \
        | assemble \
        | predict_genes

}
