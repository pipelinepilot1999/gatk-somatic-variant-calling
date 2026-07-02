#!/bin/bash
cd ~/gatk-somatic-variant-calling
REF=data/reference/Homo_sapiens_assembly38.fasta
DS=data/funcotator/funcotator_dataSources.v1.8.hg38.20230908s

gatk Funcotator \
  -R $REF \
  -V results/variants/somatic.pass.vcf.gz \
  -O results/variants/somatic.funcotated.maf \
  --output-file-format MAF \
  --data-sources-path $DS \
  --ref-version hg38
