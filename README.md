# RNA-seq Differential Expression Pipeline
### Luminal A vs Triple Negative Breast Cancer (TNBC) | TCGA-BRCA

[![Snakemake](https://img.shields.io/badge/snakemake-9.22.0-brightgreen)](https://snakemake.readthedocs.io)
[![AWS](https://img.shields.io/badge/cloud-AWS%20EC2-orange)](https://aws.amazon.com)
[![Docker](https://img.shields.io/badge/container-Docker-blue)](https://docker.com)

## Overview

End-to-end reproducible RNA-seq pipeline for differential gene expression analysis, built with Snakemake, deployed on AWS EC2, and validated on mouse mammary gland data (GSE60450) before scaling to TCGA-BRCA human breast cancer data.

**Biological Question:** What genes are differentially expressed between Luminal A and Triple Negative Breast Cancer (TNBC) subtypes?

**Clinical Relevance:** TNBC has the worst prognosis among breast cancer subtypes and lacks targeted therapies. Identifying molecular differences between Luminal A and TNBC can inform therapeutic strategies.

---

## Pipeline Overview

```
Raw FASTQ -> FastQC -> Trim Galore -> STAR -> featureCounts -> DESeq2 -> Plots
```

| Step | Tool | Purpose |
|------|------|---------|
| QC | FastQC + MultiQC | Read quality assessment |
| Trimming | Trim Galore | Adapter removal |
| Alignment | STAR 2.7.11b | Splice-aware alignment to reference genome |
| Quantification | featureCounts | Gene-level read counting |
| Differential Expression | DESeq2 | Statistical testing |
| Visualization | ggplot2 + pheatmap | Volcano, PCA, heatmap plots |

---

## Results - Phase 1 Validation (GSE60450)

Pipeline validated on mouse mammary gland dataset (Law et al. 2016) comparing basal vs luminal cell types.

### PCA Plot
![PCA Plot](results/deseq2/pca_plot.png)

PC1 captures 88% of variance with complete separation between basal and luminal samples, confirming strong transcriptional differences between cell types.

### Volcano Plot
![Volcano Plot](results/deseq2/volcano_plot.png)

Thousands of significant DEGs (padj < 0.05, |log2FC| > 1) with symmetric up and downregulation between groups.

### Heatmap - Top 50 DEGs
![Heatmap](results/deseq2/heatmap.png)

Perfect hierarchical clustering separating basal and luminal samples. Top luminal markers include known milk protein genes (Wap, Csn2, Glycam1) and basal markers include myoepithelial genes (Krt6a, Piezo2), consistent with published results from Law et al. 2016.

---

## Repository Structure

```
rnaseq-brca-pipeline/
├── workflow/
│   ├── Snakefile          # Pipeline rules
│   ├── config.yaml        # Parameters and paths
│   ├── scripts/
│   │   └── deseq2.R       # DESeq2 analysis and plots
│   └── envs/              # Conda environment files
├── data/
│   ├── raw/               # FASTQ files and metadata
│   ├── genome/            # Reference genome and STAR index
│   └── annotation/        # GTF annotation files
├── results/
│   ├── fastqc/            # FastQC HTML reports
│   ├── trimmed/           # Adapter-trimmed reads
│   ├── star/              # BAM alignment files
│   ├── counts/            # featureCounts matrix
│   ├── multiqc/           # MultiQC combined report
│   └── deseq2/            # DESeq2 results and plots
├── logs/                  # Per-rule log files
└── docker/                # Dockerfile (coming soon)
```

---

## Quickstart

### Requirements
- Snakemake 9.22.0
- AWS EC2 r5.large or larger (16GB+ RAM for STAR index)
- Conda/Mamba

### Setup

```bash
git clone git@github.com:pipelinepilot1999/rnaseq-brca-pipeline.git
cd rnaseq-brca-pipeline
nano workflow/config.yaml
conda create -c conda-forge -c bioconda -n snakemake snakemake -y
conda activate snakemake
```

### Run Pipeline

```bash
# Dry run first
snakemake -n --snakefile workflow/Snakefile

# Full run
snakemake --snakefile workflow/Snakefile --cores 8
```

---

## Project Phases

| Phase | Dataset | Contrast | Status |
|-------|---------|----------|--------|
| Phase 1 | GSE60450 (mouse mammary) | Basal vs Luminal | Complete |
| Phase 2 | TCGA-BRCA | Tumor vs Normal | In progress |
| Phase 3 | TCGA-BRCA | Luminal A vs TNBC | Planned |

---

## Tools and Technologies

- **Workflow:** Snakemake 9.22.0
- **Cloud:** AWS EC2 (r5.2xlarge for index build, t3.medium for setup)
- **Alignment:** STAR 2.7.11b (splice-aware)
- **Quantification:** featureCounts v2.1.1
- **DE Analysis:** DESeq2 (R/Bioconductor)
- **Visualization:** ggplot2, pheatmap
- **QC:** FastQC, MultiQC, Trim Galore
- **Container:** Docker (in progress)

---

## Author

**Sathvik Sai Appagana**
MS Bioinformatics, Indiana University Indianapolis
GitHub: [pipelinepilot1999](https://github.com/pipelinepilot1999)
Portfolio: [atac-seq-pipeline](https://github.com/pipelinepilot1999/atac-seq-pipeline)
