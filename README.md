# TFM
TFM Master Bioinformática y Bioestadística UOC
Febrero - Junio 2022 

Nextflow-NGS
============
[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.24.0-brightgreen.svg)](http://nextflow.io)
Flujo de trabajo para el análisis (pipeline) de datos NGS basado en [Nextflow](https://www.nextflow.io/) para el estudio de la neurofibromatosis tipo NF1.



-------------
Preprocesado 
-------------
* **RECORTE_ADAPTADORES_QC** 
* **INDICE_REFERENCIA** 
* **MAPA_A_REFERENCIA**
* **CONVERSION_ SAM_A_BAM**
* **ORDENAMIENTO**
* **CONTEO_DE_LECTURAS** 
* **MARCAR_DUPLICADOS**
* **BQSR**



----------------------------
Descubrimiento de Variantes 
----------------------------
* **DESCUBRIMIENTO_DE_VARIANTES**
* **CATEGORIZAR_VARIANTES** 


-----------------------
Anotación de Variantes 
-----------------------
* **FILTRADO_DURO**

-----------------------
ADICIONALES: 
-----------------------

* **SRA_CONVERSION** - Conversion de archivos  \*.sra del repositorio [Read Sequencing Archive](https://www.ncbi.nlm.nih.gov/sra) a formato fastq utilizando fastq-dump del kit [sra-tools](https://github.com/ncbi/sra-tools).

