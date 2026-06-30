# Germline Variant Analysis Pipeline

A comprehensive bioinformatics workflow for identifying, annotating, filtering, and prioritizing germline genetic variants from Next-Generation Sequencing (NGS) data. This project demonstrates an end-to-end variant analysis pipeline following widely accepted best practices for genomic data processing and variant interpretation.

---

## Project Overview

The objective of this project is to process raw sequencing data and identify high-confidence germline variants that may contribute to inherited genetic disorders. The workflow includes sequence quality assessment, alignment, variant calling, annotation, filtering, and clinical interpretation.

This project showcases practical experience with NGS analysis, variant interpretation, and genomic databases commonly used in research and clinical genomics.

---

## Workflow

```text
Raw FASTQ Files
        │
        ▼
Quality Control (FastQC)
        │
        ▼
Read Trimming (Fastp/Trimmomatic)
        │
        ▼
Alignment to Reference Genome (BWA-MEM)
        │
        ▼
Sorting & Indexing (SAMtools)
        │
        ▼
Duplicate Marking (Picard)
        │
        ▼
Base Quality Score Recalibration (GATK)
        │
        ▼
Variant Calling (GATK HaplotypeCaller)
        │
        ▼
Variant Filtering
        │
        ▼
Variant Annotation (ANNOVAR / VEP)
        │
        ▼
Clinical Interpretation
```

---

## Features

- End-to-end germline variant calling pipeline
- Quality control of raw sequencing reads
- High-quality read alignment
- Variant calling using GATK Best Practices
- Functional annotation of SNPs and InDels
- Population frequency filtering
- Clinical significance annotation
- Candidate variant prioritization
- Reproducible analysis workflow

---

## Tools and Software

| Tool | Purpose |
|------|----------|
| FastQC | Read quality assessment |
| Fastp / Trimmomatic | Read trimming |
| BWA-MEM | Sequence alignment |
| SAMtools | BAM processing |
| Picard | Duplicate marking |
| GATK | Variant calling & filtering |
| ANNOVAR / VEP | Variant annotation |
| bcftools | Variant manipulation |
| IGV | Variant visualization |
| MultiQC | QC report aggregation |

---

## Databases Used

- ClinVar
- dbSNP
- gnomAD
- OMIM
- RefSeq
- Ensembl
- UCSC Genome Browser

---

## Pipeline Steps

### 1. Quality Control

- Assess sequencing quality
- GC content analysis
- Adapter contamination
- Per-base sequence quality

**Output**

- FastQC reports
- MultiQC summary

---

### 2. Read Preprocessing

- Adapter removal
- Low-quality base trimming
- Filtering short reads

---

### 3. Sequence Alignment

Reads are aligned against the human reference genome using **BWA-MEM**.

Outputs:

- Sorted BAM
- Indexed BAM

---

### 4. Post Alignment Processing

- Duplicate marking
- BAM indexing
- Base quality recalibration

---

### 5. Variant Calling

Variants are identified using **GATK HaplotypeCaller**.

Generated files:

- Raw VCF
- Filtered VCF

---

### 6. Variant Filtering

Common filtering criteria include:

- Read Depth (DP)
- Quality Score (QUAL)
- Genotype Quality (GQ)
- Variant Quality Score
- Allele Balance

---

### 7. Variant Annotation

Variants are annotated using genomic databases to identify:

- Gene information
- Protein changes
- Functional consequences
- Population frequencies
- Clinical significance
- Disease associations

---

### 8. Variant Prioritization

Candidate variants are prioritized based on:

- Rare allele frequency
- Predicted pathogenicity
- ClinVar classification
- Disease relevance
- Functional impact
- Literature evidence

---

## Example Directory Structure

```text
Germline-Variant-Analysis/
│
├── data/
│   ├── raw_fastq/
│   ├── trimmed/
│   ├── bam/
│   ├── vcf/
│   └── annotation/
│
├── scripts/
│   ├── qc.sh
│   ├── alignment.sh
│   ├── variant_calling.sh
│   ├── annotation.sh
│   └── filtering.sh
│
├── results/
│   ├── fastqc/
│   ├── multiqc/
│   ├── bam/
│   ├── variants/
│   └── reports/
│
├── figures/
├── README.md
└── requirements.txt
```

---

## Results

The pipeline successfully performs:

- High-quality read processing
- Accurate alignment to the reference genome
- Germline SNP and InDel detection
- Functional variant annotation
- Identification of clinically relevant variants
- Generation of analysis-ready annotated VCF files

---

## Skills Demonstrated

- Bioinformatics
- Next-Generation Sequencing (NGS)
- Germline Variant Analysis
- Variant Annotation
- Clinical Genomics
- Bash Scripting
- Linux
- Genomic Data Processing
- Quality Control
- Variant Interpretation
- Data Analysis

---

## Future Improvements

- Structural Variant Detection
- Copy Number Variation Analysis
- Trio Analysis
- ACMG Classification Automation
- Workflow automation using Nextflow
- Docker containerization
- Cloud deployment
- Interactive reporting dashboard

---

## References

- GATK Best Practices
- BWA-MEM
- SAMtools
- ClinVar
- gnomAD
- ANNOVAR
- Ensembl Variant Effect Predictor (VEP)

---

## Author

**Abhimanyu Mandal**

Computational Biologist | Bioinformatics | NGS Analysis | Clinical Genomics | Data Analytics

LinkedIn: *(Add your profile)*

GitHub: *(Add your GitHub profile)*

---
