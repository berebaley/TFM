#!/usr/bin/env nextflow
params.readDir = "./*.sra"
params.output = "./results"
log.info " SRA CONVERSION"
log.info "========================="
log.info "readDir              : ${params.readDir}"

Channel
    .fromFilePairs( params.readDir, size: 1 )
    .ifEmpty { error "Cannot find any reads matching: ${params.readDir}" }
    .set { read_files } 

process Conversion_SRA {
    publishDir "${params.output}", mode: 'copy'
	tag { name }

	module params.sra_tk
	
	cpus = 1
	     
    input:
    set val(name), file(reads) from read_files
    
    output:
    file '*fastq.gz' into trimmed_reads
 
    script:
    """
    fastq-dump --split-3 --gzip ./$reads
    """
}
 
workflow.onComplete { 
	log.info ( workflow.success ? "Flujo de Trabajo Completado!" : "No fue posible completar el flujo de trabajo" )
}
