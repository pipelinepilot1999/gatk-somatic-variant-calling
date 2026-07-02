#!/bin/bash
cd ~/gatk-somatic-variant-calling
REF=data/reference/Homo_sapiens_assembly38.fasta
DBSNP=data/reference/Homo_sapiens_assembly38.dbsnp138.vcf
MILLS=data/reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

# Pass 1: build the recalibration table
gatk BaseRecalibrator \
  -I results/markdup/tumor.markdup.bam \
  -R $REF \
  --known-sites $DBSNP \
  --known-sites $MILLS \
  -O results/bqsr/tumor.recal.table

# Pass 2: apply the table to rewrite quality scores
gatk ApplyBQSR \
  -I results/markdup/tumor.markdup.bam \
  -R $REF \
  --bqsr-recal-file results/bqsr/tumor.recal.table \
  -O results/bqsr/tumor.recal.bam
