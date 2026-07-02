#!/bin/bash
cd ~/gatk-somatic-variant-calling
REF=data/reference/Homo_sapiens_assembly38.fasta
gatk FilterMutectCalls \
  -R $REF \
  -V results/variants/somatic.vcf.gz \
  -O results/variants/somatic.filtered.vcf.gz
