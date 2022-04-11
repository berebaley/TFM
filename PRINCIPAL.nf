
params.reads = "$baseDir/datos/reads/*_{1,2}.fq"
params.genome = "$baseDir/datos/genome.fa"
params.dbsnp = "$baseDir/datos/known_variants.vcf"
params.outdir = 'RESULTADOS'

log.info """\
         NF1-NGS   P I P E L I N E    
         =============================
         genome: ${params.genome}
         reads : ${params.reads}
         outdir: ${params.outdir}
         """
         .stripIndent()
 

/*
 * 
 */
num_samples = 0 
Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}"  }
    .tap { read_pairs_ch }
    .subscribe({ num_samples += 1 })

/*
 * 
 */
process INDICE_REFERENCIA {
    tag "$genome.baseName"
    
    input:
    path genome from params.genome
     
    output:
    path 'genome.index*' into index_ch
       
    """
    bowtie2-build --threads ${task.cpus} ${genome} genome.index
    """
}
 
/*
 * 
 */
process MAPA_A_REFERENCIA {
    tag "$pair_id"
    publishDir "${params.outdir}/bowtie_logs", pattern: '*.log'

    input:
    file index from index_ch.collect()
    set pair_id, file(reads) from read_pairs_ch


    output:
    set val(pair_id), file("${pair_id}.sam") into aligned_sam
    file "${pair_id}.log" into bowtie_logs

    script:
    """
    bowtie2 \\
        --threads ${task.cpus} \\
        -x genome.index\\
        -q -1 ${reads[0]} -2 ${reads[1]} \\
        --very-sensitive-local \\
        -S ${pair_id}.sam \\
        --no-unal \\
        2>&1 | tee ${pair_id}.log
    """

}
  
/*
 * 
 */
process CONVERSION_ SAM_A_BAM {
    tag "$pair_id"
    publishDir "${params.outdir}/bam", mode:'copy'

    input:
    set val(pair_id), file(sam) from aligned_sam

    output:
    set val(pair_id), file("${pair_id}.bam") into bam_ch, bam_stats_ch, input2_ch


    script:
    """
    samtools view -S -b $sam > ${pair_id}.bam
    """
}

/*
 * 
 */
process ORDENAMIENTO {
    tag "$pair_id"
    publishDir "${params.outdir}/samtools", mode:'copy'

    input:
    set val(pair_id), file(bam) from bam_stats_ch

    output:
    set val(pair_id),file ("${pair_id}.sorted_rg.bam") into sorted_bam
    file ("${pair_id}.sorted_rg.bam.bai") into index_bam
    file ("*stat*") into samtools_stats

    script:
    """
    samtools sort \\
            $bam \\
            -@ ${task.cpus} \\
            -o ${pair_id}.sorted.bam
    samtools addreplacerg -r ID:1 -r LB:1 -r SM:SAMPLE_${pair_id} -r PL:ILLUMINA -r PU:1 -o ${pair_id}.sorted_rg.bam ${pair_id}.sorted.bam
    samtools index ${pair_id}.sorted_rg.bam > ${pair_id}.sorted_rg.bam.bai
    samtools flagstat ${pair_id}.sorted_rg.bam > ${pair_id}.flagstats
    samtools stats ${pair_id}.sorted_rg.bam > ${pair_id}.stats
    samtools idxstats ${pair_id}.sorted_rg.bam > ${pair_id}.idxstats
    """
}

/*
 * 
 */
process CONTEO_LECTURAS {
  tag "$pair_id"
  publishDir "${params.outdir}/counts", mode:'copy'

  input:
  set val(pair_id), file(bam) from bam_ch

  output:
  file "*bowtie*" into bowtie_stats

  shell:
  '''
  samtools view -f 0x2 !{bam} | grep -v "XS:i:" | foo.py | cut -f 3 | sort | uniq -c | awk '{printf("%s\t%s\\n", $2, $1)}' > !{pair_id}_bowtie_csp_counts.txt
  samtools view -F 260 !{bam} | cut -f 3 | sort | uniq -c | awk '{printf("%s\t%s\\n", $2, $1)}'  > !{pair_id}_bowtie_all_counts.txt
  samtools view -f 0x2 !{bam} | grep -v "XS:i:" | cut -f 3 | sort | uniq -c | awk '{printf("%s\t%s\\n", $2, $1)}' > !{pair_id}_bowtie_cs_counts.txt
  '''
}

/*
 * 
 */
process MARCAR_DUPLICADOS {
    tag "$pair_id"
    publishDir "${params.outdir}/duplicates", mode:'copy'
    
    input:
    set val(pair_id), file(bam) from sorted_bam

    output:
    set val(pair_id), file ("${pair_id}.dedup.bam") into deduplicated_bams

    script:
    """
    gatk MarkDuplicates\
    --INPUT $bam\
    --OUTPUT ${pair_id}.dedup.bam\
    --METRICS_FILE ${pair_id}.metrics.txt
    """
}

/*
 * 
 */
process BQSR {
    publishDir "${params.outdir}/bqsr", mode: "copy"
    tag "${pair_id}"

    input:
    tuple val(pair_id), file(bam) from deduplicated_bams
    file(bai) from index_bam

    output:
    tuple val(pair_id), file("${pair_id}.preprocessed.bam") into recalibrated_bams
    file "${pair_id}.preprocessed.bam"
    file "${pair_id}.preprocessed.bai"

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
    """
}

/*
 * 
 */
process DESCUBRIMIENTO_VARIANTES {
    publishDir "${params.outdir}/vcf", mode: "copy"
    tag "${pair_id}"

    input:
    tuple val(pair_id), file(bam) from recalibrated_bams
    
    output:
    tuple val(pair_id), file("${pair_id}.g.vcf") into called_variants
    file "${pair_id}.g.vcf.idx"
    file "${pair_id}.out.bam"
    file "${pair_id}.out.bai"

    """ 
    gatk HaplotypeCaller -R ${params.genome} -I ${bam} -O ${pair_id}.g.vcf -bamout ${pair_id}.out.bam -ERC GVCF
    """
}

process CATEGORIZAR_VARIANTES {
    publishDir "${params.outdir}/Categorized_Variants", mode: "copy"
    tag "${pair_id}"

    input:
    tuple val(pair_id), file(vcf) from called_variants
    
    output:
    tuple val(pair_id), file("${pair_id}_snps.vcf") into snp_vcf
    tuple val(pair_id), file("${pair_id}_indels.vcf") into indels_vcf
    file "${pair_id}_snps.vcf.idx"
    file "${pair_id}_indels.vcf.idx"
    """ 
    gatk SelectVariants -V ${vcf} -select-type SNP -O ${pair_id}_snps.vcf
    gatk SelectVariants -V ${vcf} -select-type INDEL -O ${pair_id}_indels.vcf
    """
}

/*
 * 
 */
process FILTRADO_DURO {
    publishDir "${params.outdir}/Hard_Filtered", mode: "copy"
    tag "${pair_id}"

    input:
    tuple val(pair_id), file(snps) from snp_vcf
    tuple val(pair_id), file(indels) from indels_vcf
    
    output:
    tuple val(pair_id), file("${pair_id}_snps_filtered.vcf") into filtered_snp
    tuple val(pair_id), file("${pair_id}_indels_filtered.vcf") into filtered_indels
    file "${pair_id}_snps_filtered.vcf.idx"
    file "${pair_id}_indels_filtered.vcf.idx"
    
    """ 
    gatk VariantFiltration -V ${indels}  -filter "QD<2.0" --filter-name "QD2" -filter "QUAL<30.0" --filter-name "QUAL30" -filter "FS>200.0" --filter-name "FS200" -filter "ReadPosRankSum<-20.0" --filter-name "ReadPosRankSum-20" -O ${pair_id}_indels_filtered.vcf
    gatk VariantFiltration -V ${snps}  -filter "QD<2.0" --filter-name "QD2" -filter "QUAL<30.0" --filter-name "QUAL30" -filter "SOR>3.0" --filter-name "SOR3" -filter "FS>60.0" --filter-name "FS60" -filter "MQ<40.0" --filter-name "MQ40" -filter "MQRankSum<-12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum<-8.0" --filter-name "ReadPosRankSum-8" -O ${pair_id}_snps_filtered.vcf 
    """
}

workflow.onComplete { 
	log.info ( workflow.success ? "Flujo de Trabajo Completado!" : "No fue posible completar el flujo de trabajo" )
}
