params.files = "/Users/alonsobajo/berebaley/resultados_preprocesamiento/bqsr/*.bam"
params.genome = "$baseDir/hg38_v0_Homo_sapiens_assembly38.fasta"
params.sources = "$baseDir/funcotator_dataSources.v1.6.20190124g"
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

process MPILEUP {
    tag "${bam.baseName}"

    input:
    path(bam) from read_pairs_ch

    output:
    file("${bam.baseName}.pileup") into mpileup_ch

    """ 
    samtools mpileup -f ${params.genome} -q 1 -B ${bam} > ${bam.baseName}.pileup
    """
}

process VARSCAN {
    tag "${pileup.baseName}"

    input:
    file (pileup) from mpileup_ch

    output:
    file("${pileup.baseName}_output_snp_germline.vcf") into variants_snp_germline_ch
    file("${pileup.baseName}_output_indel_germline.vcf") into variants_indel_germline_ch
    file("${pileup.baseName}_output_cnc_germline.vcf") into variants_cnc_germline_ch
    file("${pileup.baseName}_output_snp_germline.vcf.idx") into variants_snp_germline_idx_ch
    file("${pileup.baseName}_output_indel_germline.vcf.idx") into variants_indel_germline_idx_ch
    file("${pileup.baseName}_output_cnc_germline.vcf.idx") into variants_cnc_germline_idx_ch
    """ 
    java -jar /Users/alonsobajo/berebaley/VarScan.v2.3.9.jar  mpileup2snp ${pileup} --output-vcf > ${pileup.baseName}_output_snp_germline.vcf
    java -jar /Users/alonsobajo/berebaley/VarScan.v2.3.9.jar  mpileup2indel ${pileup} --output-vcf > ${pileup.baseName}_output_indel_germline.vcf
    java -jar /Users/alonsobajo/berebaley/VarScan.v2.3.9.jar  mpileup2cns ${pileup} --output-vcf > ${pileup.baseName}_output_cns_germline.vcf
    gatk --java-options -Xmx16g IndexFeatureFile -I ${pileup.baseName}_output_snp_germline.vcf -O ${pileup.baseName}_output_snp_germline.vcf.idx
    gatk --java-options -Xmx16g IndexFeatureFile -I ${pileup.baseName}_output_indel_germline.vcf -O ${pileup.baseName}_output_indel_germline.vcf.idx
    gatk --java-options -Xmx16g IndexFeatureFile -I ${pileup.baseName}_output_cnc_germline.vcf -O ${pileup.baseName}_output_cnc_germline.vcf.idx

    """
}

process ANOTAR{
    publishDir "${params.outdir}/ANOTAR_GERM", mode: "copy"
    tag "${vcf1.baseName}"

    input:
    file(vcf1) from variants_snp_germline_ch
    file (ind1) from variants_snp_germline_idx_ch
    file(vcf2) from variants_indel_germline_ch
    file (ind2) from variants_indel_germline_idx_ch
    file(vcf3) from variants_cnc_germline_ch
    file (ind3) from variants_cnc_germline_idx_ch

    output:
    file "${vcf1.baseName}_snps_snpEff.ann.vcf" into snps_ann_vcf
    file "${vcf1.baseName}_snps_snpEff.ann.vcf.idx" into snps_ann_vcf_idx
    file "${vcf2.baseName}_indels_snpEff.ann.vcf" into indels_ann_vcf
    file "${vcf2.baseName}_indels_snpEff.ann.vcf.idx" into indels_ann_vcf_idx
    file "${vcf3.baseName}_cnc_snpEff.ann.vcf" into cnc_ann_vcf
    file "${vcf3.baseName}_cnc_snpEff.ann.vcf.idx" into cnc_ann_vcf_idx

    """
    java -Xmx16g -jar /Users/alonsobajo/berebaley/snpEff/snpEff.jar GRCh38.86 ${vcf1} > ${vcf1.baseName}_snps_snpEff.ann.vcf
    gatk --java-options -Xmx16g IndexFeatureFile -I ${vcf1.baseName}_snps_snpEff.ann.vcf -O ${vcf1.baseName}_snps_snpEff.ann.vcf.idx

    java -Xmx16g -jar /Users/alonsobajo/berebaley/snpEff/snpEff.jar GRCh38.86 ${vcf2} > ${vcf2.baseName}_indels_snpEff.ann.vcf
    gatk --java-options -Xmx16g IndexFeatureFile -I ${vcf2.baseName}_indels_snpEff.ann.vcf -O ${vcf2.baseName}_indels_snpEff.ann.vcf.idx

    java -Xmx16g -jar /Users/alonsobajo/berebaley/snpEff/snpEff.jar GRCh38.86 ${vcf3} > ${vcf3.baseName}_cnc_snpEff.ann.vcf
    gatk --java-options -Xmx16g IndexFeatureFile -I ${vcf3.baseName}_indels_snpEff.ann.vcf -O ${vcf3.baseName}_cnc_snpEff.ann.vcf.idx
    """

}

workflow.onComplete { 
	log.info ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
