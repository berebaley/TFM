# TFM
TFM Master Bioinformática y Bioestadística UOC
Febrero - Junio 2022 

Nextflow-NGS
============
[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.24.0-brightgreen.svg)](http://nextflow.io)
Flujo de trabajo para el análisis (pipeline) de datos NGS basado en [Nextflow](https://www.nextflow.io/) para el estudio de la neurofibromatosis tipo NF1.




Preprocesado 
-------------
* **Adapter trimming** - Standard adapter trimming and QC using [trim_galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
Genome tools
------------
* **star_align** - [VBCF](http://www.vbcf.ac.at/facilities/next-generation-sequencing/) RNA-seq alignment and QA pipeline using [STAR](https://github.com/alexdobin/STAR)
* **ChIPSeq_align** - [VBCF](http://www.vbcf.ac.at/facilities/next-generation-sequencing/) ChIP-seq alignment and QA pipeline using [bowtie](http://bowtie-bio.sourceforge.net/index.shtml)


Descubrimiento de Variantes 
---------------
* **GATK RNA-Seq** - Variant calling on RNA-Seq data using Broad's Genome Analysis Toolkit following their [best practices](https://software.broadinstitute.org/gatk/documentation/article.php?id=3891).
Peak calling
------------
* **SICER** - Standard broad peak calling using [SICER](http://home.gwu.edu/~wpeng/Software.htm).

Anotación de Variantes 
---------------


ADICIONALES: 
-----------------------
Conversion de Archivos 
-----------------------
* **SRA conversion** - Conversion de archivos  \*.sra del repositorio [Read Sequencing Archive](https://www.ncbi.nlm.nih.gov/sra) a formato fastq utilizando fastq-dump del kit [sra-tools](https://github.com/ncbi/sra-tools).

