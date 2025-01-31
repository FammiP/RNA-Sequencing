#!/bin/bash

#================================================================
# SCRIPT NAME:
#   01_FastQC.sh
#
# DESCRIPTION:
#   This script performs quality control analysis on RNA-seq reads
#   using FastQC, producing detailed reports on read quality,
#   adapter content, and other metrics. It is designed to handle
#   multiple FASTQ files in a given directory and output results
#   to a specified directory.
#
# USAGE:
#   sbatch 01_FastQC.sh
#
# REQUIREMENTS:
#   - FastQC version 0.11.9 (apptainer container specified)
#   - Access to a compute cluster with SLURM scheduler
#
# INPUT:
#   - FASTQ files containing RNA-seq reads (*.fastq.gz) located in:
#     $WORKDIR/toxoplasma_de/reads/
#
# OUTPUT:
#   - FastQC reports (HTML and zip files) in the directory:
#     $WORKDIR/FastQC
#
#================================================================

#SBATCH --job-name=FastQC
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/fparokkaran/RNA_seq/FastQC/logfile/%x_%j.out
#SBATCH --error=/data/users/fparokkaran/RNA_seq/FastQC/logfile/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=fammi.parokkaran@students.unibe.ch


WORKDIR=/data/users/fparokkaran/RNA_seq

mkdir -p $WORKDIR/FastQC

apptainer exec \
--bind $WORKDIR \
/containers/apptainer/fastqc-0.12.1.sif \
fastqc -t 8 -o $WORKDIR/O1_Quality_Control/FastQC $WORKDIR/toxoplasma_de/reads/*.fastq.gz


