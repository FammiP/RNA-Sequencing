#!/bin/bash

#================================================================
# SCRIPT NAME:
#   02_multiQC.sh
#
# DESCRIPTION:
#   This script aggregates quality control reports using MultiQC.
#   It scans the specified directory for FastQC reports and
#   compiles them into a comprehensive summary report for easy
#   assessment of the sequencing data quality.
#
# USAGE:
#   sbatch 02_multiQC.sh
#
# REQUIREMENTS:
#   - MultiQC version 1.11
#   - Existing FastQC output directory containing reports
#
# INPUT:
#   - Directory containing FastQC reports: $WORKDIR/01_Quality_Control/FastQC
#
# OUTPUT:
#   - MultiQC summary report (HTML and supporting files) in:
#     $WORKDIR/01_Quality_Control/MultiQC
#
# NOTES:
#   - Ensure FastQC has been run and reports are present in the specified directory.
#   - Update paths and email address if necessary.
#================================================================


#SBATCH --job-name=MultiQC
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/fparokkaran/RNA_seq/MultiQC/logfile/%x_%j.out
#SBATCH --error=/data/users/fparokkaran/RNA_seq/MultiQC/logfile/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=fammi.parokkaran@students.unibe.ch

ml MultiQC/1.11-foss-2021a

# Define the working directory
WORKDIR=/data/users/fparokkaran/RNA_seq

# Run MultiQC to aggregate existing FastQC reports
multiqc $WORKDIR/01_Quality_Control/FastQC -o $WORKDIR/01_Quality_Control/MultiQC