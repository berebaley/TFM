params.files = "/Users/alonsobajo/berebaley/resultados_preprocesamiento/bqsr/*.bam"

params.genome = "$baseDir/hg38_v0_Homo_sapiens_assembly38.fasta"
params.dbsnp = "$baseDir/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
params.hapmap="$baseDir/hapmap_3.3.hg38.vcf.gz"
params.omni="$baseDir/1000G_omni2.5.hg38.vcf.gz"
params.thousand="$baseDir/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
params.mills="$baseDir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

params.outdir = 'resultados_germline'

log.info """\
         NF1 NGS P I P E L I N E    
         =============================
         genome: ${params.genome}
         outdir: ${params.outdir}
         """
         .stripIndent()
 
Channel
    .fromPath(params.files)
    .ifEmpty { error "Cannot find any reads matching: ${params.files}"  }
    .set { read_pairs_ch }

process HAPLOTYPE {
    publishDir "${params.outdir}/vcf", mode: "copy"
    tag "${bam.baseName}"

    input:
    path(bam) from read_pairs_ch

    output:
    file("${bam.baseName}.g.vcf") into called_variants
    file "${bam.baseName}.g.vcf.idx" into called_variants_idx
    file "${bam.baseName}.out.bam"
    file "${bam.baseName}.out.bai"

    """ 
    gatk --java-options -Xmx16g HaplotypeCaller -R ${params.genome} -I ${bam} -O ${bam.baseName}.g.vcf -bamout ${bam.baseName}.out.bam -ERC GVCF
    """
}

process CON_GVCF{
    publishDir "${params.outdir}/db", mode: "copy"
    tag "${gvcf.baseName}"

    input:
    file (gvcf) from called_variants
    file (vcf_idx) from called_variants_idx 

    output: 
    file "${gvcf.baseName}_consolidated.vcf" into consolidated_vcf 
    file "${gvcf.baseName}_consolidated.vcf.idx" into consolidated_vcf_idx

    script:

    """
    gatk --java-options -Xmx16g GenomicsDBImport -R ${params.genome} -V ${gvcf} --genomicsdb-workspace-path my_database --intervals chr1 --intervals chr2 --intervals chr3 --intervals chr4 --intervals chr5 --intervals chr6 --intervals chr7 --intervals chr8 --intervals chr9 --intervals chr10 --intervals chr11 --intervals chr12 --intervals chr13 --intervals chr14 --intervals chr15 --intervals chr16 --intervals chr17 --intervals chr18 --intervals chr19 --intervals chr20 --intervals chr21 --intervals chr22
    gatk --java-options -Xmx16g GenotypeGVCFs -R ${params.genome} -V gendb://my_database -O ${gvcf.baseName}_consolidated.vcf
    """
}

process CATEGORIZAR_VARIANTES {
    publishDir "${params.outdir}/Categorized_Variants", mode: "copy"
    tag "${vcf.baseName}"

    input:
    file(vcf) from consolidated_vcf
    file (index) from consolidated_vcf_idx
    
    output:
    file("${vcf.baseName}_snps.vcf") into snp_vcf
    file("${vcf.baseName}_indels.vcf") into indels_vcf
    file "${vcf.baseName}_snps.vcf.idx" into snp_vcf_idx
    file "${vcf.baseName}_indels.vcf.idx" into indels_vcf_idx
    """ 
    gatk SelectVariants -V ${vcf} -select-type SNP -O ${vcf.baseName}_snps.vcf
    gatk SelectVariants -V ${vcf} -select-type INDEL -O ${vcf.baseName}_indels.vcf
    """
}

process VSQR {
    publishDir "${params.outdir}/VSQR", mode: "copy"
    tag "${snp_vcf.baseName}"

    input:
    file (snp_vcf) from snp_vcf
    file (indel_vcf) from indels_vcf
    file (snp_vcf_idx) from snp_vcf_idx
    file (indel_vcf_idx) from indels_vcf_idx
    
    output:
    file "${snp_vcf.baseName}_consolidated_recal_snps.vcf" into recal_vcf
    file "${snp_vcf.baseName}_consolidated_recal_snps.vcf.idx" into recal_vcf_idx
    file "${indel_vcf.baseName}_consolidated_recal_indels.vcf" into recal_indel_vcf
    file "${indel_vcf.baseName}_consolidated_recal_indels.vcf.idx" into recal_indel_vcf_idx

    """
    gatk --java-options -Xmx16g VariantRecalibrator -V ${snp_vcf}  --trust-all-polymorphic -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP -mode SNP --max-gaussians 6 --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${params.hapmap} --resource:omni,known=false,training=true,truth=true,prior=12.0 ${params.omni} --resource:1000G,known=false,training=true,truth=true,prior=10.0 ${params.thousand} --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.dbsnp} -O ${snp_vcf.baseName}_consolidated_snps.recal --tranches-file ${snp_vcf.baseName}_consolidated_recal_snps.tranches 
    gatk --java-options -Xmx16g ApplyVQSR -V ${snp_vcf} -O ${snp_vcf.baseName}_consolidated_recal_snps.vcf --truth-sensitivity-filter-level 99.7 --tranches-file ${snp_vcf.baseName}_consolidated_recal_snps.tranches --recal-file ${snp_vcf.baseName}_consolidated_snps.recal -mode SNP --create-output-variant-index true
    
    gatk --java-options -Xms16G VariantRecalibrator \
  -tranche 100.0 -tranche 99.95 -tranche 99.9 \
  -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
  -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 \
  -tranche 92.0 -tranche 91.0 -tranche 90.0 \
  -R ${params.genome}\
  -V ${indel_vcf} \
  --resource:mills,known=false,training=true,truth=true,prior=12.0 ${params.mills}\
  --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.dbsnp}\
  -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
  -mode INDEL -O ${indel_vcf.baseName}_merged_indel1.recal --tranches-file ${indel_vcf.baseName}_output_indel1.tranches

    gatk --java-options -Xmx16g ApplyVQSR -V ${indel_vcf} -O ${indel_vcf.baseName}_consolidated_recal_indels.vcf --truth-sensitivity-filter-level 99.7 --tranches-file ${indel_vcf.baseName}_output_indel1.tranches --recal-file ${indel_vcf.baseName}_merged_indel1.recal -mode INDEL --create-output-variant-index true
    """
}

process ANOTAR{
    publishDir "${params.outdir}/ANOTAR", mode: "copy"
    tag "${vcf1.baseName}"

    input:
    file(vcf1) from recal_vcf
    file (ind1) from recal_vcf_idx
    file(vcf2) from recal_indel_vcf
    file (ind2) from recal_indel_vcf_idx

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
