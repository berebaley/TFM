params.files = "/Users/alonsobajo/berebaley/resultados_preprocesamiento/bqsr/*_{1,2}.bam"
params.genome = "$baseDir/hg38_v0_Homo_sapiens_assembly38.fasta"
params.germline = "$baseDir/af-only-gnomad.hg38.vcf.gz"
params.normals = "$baseDir/1000g_pon.hg38.vcf.gz"
params.small = "$baseDir/small_exac_common_3.hg38.vcf.gz"
params.sources = "$baseDir/funcotator_dataSources.v1.6.20190124g"
params.outdir = 'resultados_somatic'

log.info """\
         NF1 NGS P I P E L I N E    
         =============================
         genome: ${params.genome}
         outdir: ${params.outdir}
         """
         .stripIndent()
 
Channel
    .fromFilePairs(params.files)
    .ifEmpty { error "Cannot find any reads matching: ${params.files}"  }
    .set { read_pairs_ch}

Channel
    .fromFilePairs(params.files)
    .ifEmpty { error "Cannot find any reads matching: ${params.files}"  }
    .set { read_pairs_ch2}

process MUTECT {
    publishDir "${params.outdir}/Mutect", mode: "copy"
    tag "${pair_id}"

    input:
     tuple val(pair_id), path(reads) from read_pairs_ch

    output:
    tuple val(pair_id), file("${pair_id}_somatic.vcf.gz") into called_variants
    file "${pair_id}_somatic.vcf.gz.tbi" into called_variants_idx
    file "${pair_id}_somatic.vcf.gz.stats" into called_variants_stats

    """ 
    gatk --java-options -Xmx16g Mutect2 -R ${params.genome} -I ${reads[0]} -tumor SAMPLE_${pair_id}-T-WEX.read -I ${reads[1]} -normal SAMPLE_${pair_id}-N-WEX.read --germline-resource ${params.germline} --panel-of-normals ${params.normals} -O ${pair_id}_somatic.vcf.gz
    """
}

process GETPILEUPSUM {
    publishDir "${params.outdir}/gpus", mode: "copy"
    tag "${pair_id}"

    input:
    tuple val(pair_id), path(reads) from read_pairs_ch2

    output:
    tuple val(pair_id), file("${pair_id}_tumor_getpileupsummaries.table"),  file("${pair_id}_normal_getpileupsummaries.table") into gpus_ch

    """ 
    gatk GetPileupSummaries -I ${reads[0]} -V ${params.small} -L ${params.small} -O ${pair_id}_tumor_getpileupsummaries.table
    gatk GetPileupSummaries -I ${reads[1]} -V ${params.small} -L ${params.small} -O ${pair_id}_normal_getpileupsummaries.table
    """
}
process CALC_CONTAMINATION {
    publishDir "${params.outdir}/contamination", mode: "copy"
    tag "${pair_id}"

    input:
    tuple val(pair_id), file (tab1) , file (tab2) from gpus_ch 

    output:
    tuple val(pair_id), file("${pair_id}_contamination.table") into contamination_ch

    """ 
    gatk CalculateContamination  -I ${tab1} -matched ${tab2} -O ${pair_id}_contamination.table                   
    """
}

process FILTER_SOMATIC {
    publishDir "${params.outdir}/filter_somatic", mode: "copy"
    tag "${pair_id}"

    input:
    tuple val(pair_id), file(vcf) from called_variants
    file (vcf_idx) from called_variants_idx
    tuple val(pair_id), file(table) from contamination_ch
    file (stats) from called_variants_stats

    output:
    tuple val(pair_id), file("${pair_id}_somatic_filtered.vcf.gz") into filter_somatic_ch
    file "${pair_id}_somatic_filtered.vcf.gz.tbi" into filter_somatic_idx

    """ 
    gatk FilterMutectCalls -R ${params.genome} -V ${vcf} --contamination-table ${table} --stats ${stats} -O ${pair_id}_somatic_filtered.vcf.gz                 
    """
}


process FUNCOTATOR {
    publishDir "${params.outdir}/funcotator", mode: "copy"
    tag "${pair_id}"

    input:
    tuple val(pair_id), file (vcf) from filter_somatic_ch
    file (vcf_idx) from filter_somatic_idx

    output:
    tuple val(pair_id), file("${pair_id}_funcotate_somatic.vcf.gz") into annotated_somatic_ch
    file "${pair_id}_funcotate_somatic.vcf.gz.tbi" into annotated_somatic_idx

    """ 
    gatk Funcotator --data-sources-path ${params.sources} --ref-version hg38 -R ${params.genome} -V ${vcf} -O ${pair_id}_funcotate_somatic.vcf.gz --output-file-format VCF                   
    """
}

workflow.onComplete { 
	log.info ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
