#!/usr/bin/env nextflow
params.memory = "3g"
params.cpus = 1
params.readDir = "./results/SRR2079547.fastq.gz"
params.output = "./resultsfastp"
log.info " FASTP pipeline "
log.info "========================="
log.info "readDir              : ${params.readDir}"

Channel
    .fromFilePairs( params.readDir, size: 1 )
    .ifEmpty { error "Cannot find any reads matching: ${params.readDir}" }
    .set { fastq1 } 


process RECORTE_ADAPTADORES_QC {

    cpus params.cpus
    memory params.memory
    publishDir "${params.output}", mode: "copy"
    tag "${name}"

    conda (params.enable_conda ? "bioconda::fastp=0.20.1" : null)

    input:
        tuple val(name), file(fastq_reads) from fastq1

    output:
        tuple val(name), file("${fastq_reads.baseName}.trimmed.fq.gz")
        file("${name}.fastp_stats.json")
        file("${name}.fastp_stats.html")

    """
    # --input_files needs to be forced, otherwise it is inherited from profile in tests
    fastp \
    --in1 ${fastq_reads} \
    --out1 ${fastq_reads.baseName}.trimmed.fq.gz \
    --json ${name}.fastp_stats.json \
    --html ${name}.fastp_stats.html
    """
}