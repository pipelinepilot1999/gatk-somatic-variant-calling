#!/bin/bash
cd ~/gatk-somatic-variant-calling
gatk MarkDuplicates \
  -I results/aligned/normal.sorted.bam \
  -O results/markdup/normal.markdup.bam \
  -M results/markdup/normal.markdup.metrics.txt
