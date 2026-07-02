#!/bin/bash
cd ~/gatk-somatic-variant-calling
REF=data/reference/Homo_sapiens_assembly38.fasta
bwa-mem2 mem -t 16 \
  -R '@RG\tID:normal\tSM:normal\tPL:ILLUMINA\tLB:normal' \
  $REF \
  data/raw/tcrb/fastq/SRR2089359_1.fastq \
  data/raw/tcrb/fastq/SRR2089359_2.fastq \
  | samtools sort -@ 4 -o results/aligned/normal.sorted.bam -
