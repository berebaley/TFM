params.files = "/Users/alonsobajo/berebaley/resultados_preprocesamiento/bqsr/*_{1,2}.bam"
params.genome = "$baseDir/hg38_v0_Homo_sapiens_assembly38.fasta"
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


process MPILEUP {
    tag "${pair_id}"

    input:
    tuple val(pair_id), path(reads) from read_pairs_ch

    output:
    tuple val(pair_id), file("${pair_id}.pileup") into mpileup_ch

    """ 
    samtools mpileup -f ${params.genome} -q 1 -B ${reads[1]} ${reads[0]} > ${pair_id}.pileup
    """
}

process VARSCAN {
    tag "${pair_id}"

    input:
    tuple val(pair_id), file (pileup) from mpileup_ch

    output:
    tuple val(pair_id), file("${pair_id}_output.snp.vcf"), file("${pair_id}_output.indel.vcf") into variants_somatic_ch

    """ 
    java -jar /Users/alonsobajo/berebaley/VarScan.v2.3.9.jar somatic ${pileup} ${pair_id}_output --mpileup 1 --min-coverage 8 --min-coverage-normal 8 --min-coverage-tumor 6 --min-var-freq 0.10 --min-freq-for-hom 0.75 --normal-purity 1.0 --tumor-purity 1.00 --p-value 0.99 --somatic-p-value 0.05 --strand-filter 0 --output-vcf
    """
}

process PROCESAR_VARIANTES {
    publishDir "${params.outdir}/process_var", mode: "copy"
    tag "${pair_id}"

    input:
    tuple val(pair_id), file (var1), file (var2) from variants_somatic_ch
    
    output:
    tuple val(pair_id), file("${pair_id}_output.snp.Somatic.vcf") into process_variants_somatic_snp_ch
    tuple val(pair_id), file("${pair_id}_output.indel.Somatic.vcf") into process_variants_somatic_indel_ch
    file "${pair_id}_output.snp.Somatic.vcf.idx" into process_variants_somatic_snp_idx
    file "${pair_id}_output.indel.Somatic.vcf.idx" into process_variants_somatic_indel_idx

    """ 
    java -jar /Users/alonsobajo/berebaley/VarScan.v2.3.9.jar processSomatic ${var1} --min-tumor-freq 0.10 --max-normal-freq 0.05 --p-value 0.07
    java -jar /Users/alonsobajo/berebaley/VarScan.v2.3.9.jar processSomatic ${var2} --min-tumor-freq 0.10 --max-normal-freq 0.05 --p-value 0.07
    gatk --java-options -Xmx16g IndexFeatureFile -I ${pair_id}_output.snp.Somatic.vcf -O ${pair_id}_output.snp.Somatic.vcf.idx
    gatk --java-options -Xmx16g IndexFeatureFile -I ${pair_id}_output.indel.Somatic.vcf -O ${pair_id}_output.indel.Somatic.vcf.idx
    
    """
}
process ANOTAR{
    publishDir "${params.outdir}/ANOTAR", mode: "copy"
    tag "${pair_id}"

    input:
    tuple val(pair_id),file(vcf1) from process_variants_somatic_snp_ch
    file (ind1) from process_variants_somatic_snp_idx
    tuple val(pair_id),file(vcf2) from process_variants_somatic_indel_ch
    file (ind2) from process_variants_somatic_indel_idx

    output:
    file "${vcf1.baseName}_snps_snpEff.ann.vcf" into snps_ann_vcf
    file "${vcf1.baseName}_snps_snpEff.ann.vcf.idx" into snps_ann_vcf_idx
    file "${vcf2.baseName}_indels_snpEff.ann.vcf" into indels_ann_vcf
    file "${vcf2.baseName}_indels_snpEff.ann.vcf.idx" into indels_ann_vcf_idx

    """
    java -Xmx16g -jar /Users/alonsobajo/berebaley/snpEff/snpEff.jar GRCh38.86 ${vcf1} > ${vcf1.baseName}_snps_snpEff.ann.vcf
    gatk --java-options -Xmx16g IndexFeatureFile -I ${vcf1.baseName}_snps_snpEff.ann.vcf -O ${vcf1.baseName}_snps_snpEff.ann.vcf.idx

    java -Xmx16g -jar /Users/alonsobajo/berebaley/snpEff/snpEff.jar GRCh38.86 ${vcf2} > ${vcf2.baseName}_indels_snpEff.ann.vcf
    gatk --java-options -Xmx16g IndexFeatureFile -I ${vcf2.baseName}_indels_snpEff.ann.vcf -O ${vcf2.baseName}_indels_snpEff.ann.vcf.idx
    """

}
workflow.onComplete { 
	log.info ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
