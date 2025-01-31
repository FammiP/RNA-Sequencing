#!/bin/bash

#================================================================
# SCRIPT NAME:
#   04_mapping2.sh
#
# DESCRIPTION:
#   This script maps paired-end RNA-seq reads to a reference genome 
#   using HISAT2. It processes multiple samples and generates SAM files 
#   for each sample. Additionally, a mapping summary log is produced 
#   for each sample.
#
# USAGE:
#   sbatch 04_mapping2.sh
#
# REQUIREMENTS:
#   - HISAT2 version 2.2.1 
#   - Indexed reference genome created with HISAT2
#   - Paired-end RNA-seq FASTQ files (*.fastq.gz)
#
# INPUT:
#   - Indexed reference genome: $REF_GEN
#   - Paired-end reads located in: $READSDIR
#   - List of sample IDs: Defined in the SAMPLES array
#
# OUTPUT:
#   - SAM files for each sample in: $OUTDIR
#   - Mapping summary logs for each sample
#
# NOTES:
#   - Ensure the reference genome is properly indexed before running the script.
#   - Update sample names, paths, and email address if necessary.
#================================================================


#SBATCH --job-name=RNASeq_Mapping
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/fparokkaran/RNA_seq/logfile/%x_%j.out
#SBATCH --error=/data/users/fparokkaran/RNA_seq/logfile/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=fammi.parokkaran@students.unibe.ch

# Define paths
WORKDIR="/data/users/fparokkaran/RNA_seq"
REF_GEN="$WORKDIR/Reference_Genome_Data/Mus_musculus"
READSDIR="$WORKDIR/Mouse_reads_data/toxoplasma_de/reads"  
OUTDIR="$WORKDIR/02_Mapping/SAM_files"


# Load Hisat2 module
ml HISAT2/2.2.1-gompi-2021a

# Create output directory if it doesn't exist
mkdir -p $OUTDIR

SAMPLES=(
    SRR7821918 SRR7821919 SRR7821920 SRR7821921 SRR7821922
    SRR7821937 SRR7821938 SRR7821939 SRR7821949 SRR7821950
    SRR7821951 SRR7821952 SRR7821953 SRR7821968 SRR7821969
    SRR7821970
)
# Loop through samples and map reads
for sample in "${SAMPLES[@]}"; do
    hisat2 -p 4 -x "$REF_GEN" \
        -1 "$READSDIR/${sample}_1.fastq.gz" \
        -2 "$READSDIR/${sample}_2.fastq.gz" \
        --rna-strandness RF \
        -S "$OUTDIR/${sample}.sam"\
        2> "$OUTDIR/${sample}_summary.log"

        
done

