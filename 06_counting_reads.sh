#!/bin/bash

#================================================================
# SCRIPT NAME:
#   06_counting_reads.sh
#
# DESCRIPTION:
#   This script uses featureCounts to count reads mapped to genomic 
#   features, such as genes, based on an annotation file. It processes 
#   sorted BAM files, generating a read count matrix for downstream 
#   differential expression analysis.
#
# USAGE:
#   sbatch 06_counting_reads.sh
#
# REQUIREMENTS:
#   - featureCounts available via the Subread module
#   - Sorted BAM files produced by previous mapping and sorting steps
#   - Annotation file in GTF format
#
# INPUT:
#   - Sorted BAM files located in: $BAMDIR
#   - Annotation file: $GTF_FILE (in GTF format)
#
# OUTPUT:
#   - Read count matrix: $OUTDIR/fCounts_any_end_aligned.txt
#
# NOTES:
#   - Ensure the BAM files are properly sorted and indexed before running the script.
#   - Update paths, node allocation, and email address if necessary.
#================================================================


#SBATCH --job-name=Read_Counting
#SBATCH --time=13:15:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --partition=pibu_el8
#SBATCH --nodelist=binfservas16
#SBATCH --exclude=binfservas32
#SBATCH --output=/data/users/fparokkaran/RNA_seq/logfile/%x_%j.out
#SBATCH --error=/data/users/fparokkaran/RNA_seq/logfile/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=fammi.parokkaran@students.unibe.ch

# Define paths
WORKDIR="/data/users/fparokkaran/RNA_seq"
BAMDIR="$WORKDIR/02_Mapping/BAM_Files"      # Directory containing sorted BAM files
OUTDIR="$WORKDIR/03_FeatureCounts"         # Directory to store read count results
GTF_FILE="$WORKDIR/Reference_Genome_Data/annotation/Mus_musculus.GRCm39.113.gtf"  # Path to the annotation file

# Load Subread module (featureCounts is part of Subread)
ml Subread/2.0.3-GCC-10.3.0

# Create output directory if it doesn't exist
mkdir -p $OUTDIR

# Perform read counting with featureCounts
featureCounts -T 16 \
    -a "$GTF_FILE" \
    -o "$OUTDIR/fCounts_any_end_aligned.txt" \
    -p \
    --countReadPairs \
    "$BAMDIR"/*_sorted.bam

