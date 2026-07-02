#!/bin/bash
cd ~/gatk-somatic-variant-calling
gatk MarkDuplicates \
  -I results/aligned/tumor.sorted.bam \
  -O results/markdup/tumor.markdup.bam \
  -M results/markdup/tumor.markdup.metrics.txt
