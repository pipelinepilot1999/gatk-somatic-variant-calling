#!/bin/bash
cd ~/gatk-somatic-variant-calling
REF=data/reference/Homo_sapiens_assembly38.fasta
gatk HaplotypeCaller \
  -I results/bqsr/normal.recal.bam \
  -R $REF \
  -O results/variants/normal.germline.vcf.gz
