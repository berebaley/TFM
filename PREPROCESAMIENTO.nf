params.reads = "$baseDir/SRA/*_{1,2}.fastq.gz"
params.genome = "$baseDir/hg38_v0_Homo_sapiens_assembly38.fasta"
params.outdir = 'resultados_preprocesamiento'
params.dbsnp = "$baseDir/Homo_sapiens_assembly38.dbsnp138.vcf.gz"

log.info """\
         NF1 NGS P I P E L I N E    
         =============================
         genome: ${params.genome}
         reads : ${params.reads}
         outdir: ${params.outdir}
         """
         .stripIndent()
 
Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}"  }
    .set { read_pairs_ch }

process RECORTE_ADAPTADORES_QC {
    
    tag "${pair_id}"

    input:
        tuple val(pair_id), path(reads) from read_pairs_ch

    output:
        tuple val(pair_id), file("${pair_id}_1.trimmed.fq.gz"), file("${pair_id}_2.trimmed.fq.gz") into trimmed_ch1
        file("${pair_id}.fastp_stats.html")

    """
    # --input_files needs to be forced, otherwise it is inherited from profile in tests
    fastp \
    -i ${reads[0]} \
    -I ${reads[1]} \
    -o ${pair_id}_1.trimmed.fq.gz \
    -O ${pair_id}_2.trimmed.fq.gz \
    --json ${pair_id}.fastp_stats.json \
    --html ${pair_id}.fastp_stats.html
    """
}


process MAPA_A_REFERENCIA {
    
    tag "${pair_id}"

    input:
        tuple val(pair_id), file (fq1), file (fq2) from trimmed_ch1

    output:
        tuple val(pair_id), file("${pair_id}.bam") into aligned_bam
        file ("${pair_id}.bam.bai") into aligned_bam_index
    script:
    """
    bwa mem -t 1 ${params.genome} ${fq1} ${fq2} | \
    samtools view -uS - | \
    samtools sort - > ${pair_id}.bam
    samtools index ${pair_id}.bam > ${pair_id}.bam.bai
    """
}

process ORDENAMIENTO {
    tag "$pair_id"
    

    input:
    set val(pair_id), file(bam) from aligned_bam
    file (idx) from aligned_bam_index
    
    output:
    set val(pair_id),file ("${pair_id}.sorted_rg.bam") into sorted_bam
    file ("${pair_id}.sorted_rg.bam.bai") into index_bam
    

    script:
    """
    samtools addreplacerg -r ID:1 -r LB:1 -r SM:SAMPLE_${pair_id} -r PL:ILLUMINA -r PU:1 -o ${pair_id}.sorted_rg.bam ${pair_id}.bam
    samtools index ${pair_id}.sorted_rg.bam > ${pair_id}.sorted_rg.bam.bai

    """
}

process MARCAR_DUPLICADOS {
    tag "$pair_id"
    
    
    input:
    set val(pair_id), file(bam) from sorted_bam
    file (idx) from index_bam

    output:
    set val(pair_id), file ("${pair_id}.dedup.bam") into deduplicated_bams
    file ("${pair_id}.dedup.bam.bai") into dedup_index_bam

    script:
    """
    gatk MarkDuplicates\
    --INPUT $bam\
    --OUTPUT ${pair_id}.dedup.bam\
    --METRICS_FILE ${pair_id}.metrics.txt 

    samtools index ${pair_id}.dedup.bam > ${pair_id}.dedup.bam.bai
    """
}

process BQSR {
    publishDir "${params.outdir}/bqsr", mode: "copy"
    tag "${pair_id}"

    input:
    tuple val(pair_id), file(bam) from deduplicated_bams
    file (idx) from dedup_index_bam

    output:
    tuple val(pair_id), file("${pair_id}.preprocessed.bam") into recalibrated_bams
    file "${pair_id}.preprocessed.bai" into recal_bam_idx

    """ 
    gatk BaseRecalibrator \
    --input ${bam} \
    --reference ${params.genome} \
    --known-sites ${params.dbsnp} \
    --output ${pair_id}.table 

    gatk ApplyBQSR \
    -R ${params.genome} \
    -I ${bam} \
    --bqsr-recal-file ${pair_id}.table \
    -O ${pair_id}.preprocessed.bam

    samtools index ${pair_id}.preprocessed.bam > ${pair_id}.preprocessed.bam.bai
    """
}

workflow.onComplete { 
	log.info ( workflow.success ? "Listo!" : "Algo salio mal" )
}
