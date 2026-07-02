#!/bin/bash
cd ~/gatk-somatic-variant-calling
REF=data/reference/Homo_sapiens_assembly38.fasta
bwa-mem2 mem -t 16 \
  -R '@RG\tID:tumor\tSM:tumor\tPL:ILLUMINA\tLB:tumor' \
  $REF \
  data/raw/tcrb/fastq/SRR2089355_1.fastq \
  data/raw/tcrb/fastq/SRR2089355_2.fastq \
  | samtools sort -@ 4 -o results/aligned/tumor.sorted.bam -
