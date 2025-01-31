# RNA-seq Analysis

This repository outlines the workflow for detecting differentially expressed genes from bulk RNA-seq data of mice infected by Toxoplasma gondii. 
The analysis workflow is part of the RNA-seq course from the University of Bern (467713).

## Installation

1. Clone this repository to your local machine.

```shell
https://github.com/FammiP/Genome-Transcriptome-Assembly-and-Annotation
```

## Step-by-Step Workflow

### Dataset
The dataset contains paired-end sequencing reads from Mus musculus lung and blood samples infected with Toxoplasma gondii. A total of 16 samples were analyzed, with biological replicates for each condition, including both infected and uninfected control samples. The samples were prepared using an Illumina platform. 

### Quality Control

#### FastQC
Run the following script to assess the quality of the raw reads
```bash
sbatch 01_FastQC.sh
```
#### MultiQC
Generates a grouped quality report using MultiQC
```bash
sbatch 02_multiQC.sh
```


### Mapping to the Reference Genome

#### Indexing 
indexes a reference genome using HISAT2 Reference genome: GRCm39 (from Ensembl)
```bash
sbatch 03_indexing_ref_gen.sh
```
#### Mapping
maps the reads to a reference genome 

```bash
sbatch 04_Mapping.sh
```

### Processing BAM Files
converts SAM files generated from RNA-seq mapping into compressed and sorted BAM files using Samtools 
```bash
sbatch 05_process_bam_files.sh
```

### Counting Reads 
Uses featureCounts to count mapped reads  
```bash
sbatch 06_counting_reads.sh
```

### Differential Expression Analysis
This script performs differential expression analysis and gene ontology (GO) enrichment for RNA-seq data using DESeq2 and clusterProfiler. It includes exploratory data visualization (PCA plots, volcano plots) and GO term visualizations (dot plots, bar plots, and tree plots)

```bash
srun --time=01:00:00 --mem=4G --cpus-per-task=1 --pty /bin/bash
Rscript DESeq2.R
```