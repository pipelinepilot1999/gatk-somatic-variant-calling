#!/bin/bash
cd ~/gatk-somatic-variant-calling
REF=data/reference/Homo_sapiens_assembly38.fasta
GNOMAD=data/reference/af-only-gnomad.hg38.vcf.gz
PON=data/reference/1000g_pon.hg38.vcf.gz

gatk --java-options "-Xmx24g" Mutect2 \
  -R $REF \
  -I results/bqsr/tumor.recal.bam \
  -I results/bqsr/normal.recal.bam \
  -normal normal \
  --germline-resource $GNOMAD \
  --panel-of-normals $PON \
  -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY  \
  -O results/variants/somatic.vcf.gz
