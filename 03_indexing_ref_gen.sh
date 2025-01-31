#!/bin/bash

#================================================================
# SCRIPT NAME:
#   03_indexing_ref_gen.sh
#
# DESCRIPTION:
#   This script indexes a reference genome using HISAT2, which is 
#   essential for efficient alignment of RNA-seq reads. The indexed 
#   files are saved in the reference genome directory for downstream 
#   processing.
#
# USAGE:
#   sbatch 03_indexing_ref_gen.sh
#
# REQUIREMENTS:
#   - HISAT2 available within the specified Apptainer container
#   - Reference genome in FASTA format
#
# INPUT:
#   - Reference genome file: Mus_musculus.GRCm39.dna.primary_assembly.fa 
#     located in $REFDIR
#
# OUTPUT:
#   - HISAT2 index files with the prefix: 
#     Mus_musculus.GRCm39.dna.primary_assembly.index
#
# NOTES:
#   - Ensure the reference genome is correctly located in $REFDIR.
#   - Update paths and email address if necessary.
#================================================================

#SBATCH --job-name=Indexing_r

#SBATCH --job-name=Indexing_reference_genome
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/fparokkaran/RNA_seq/logfile/%x_%j.out
#SBATCH --error=/data/users/fparokkaran/RNA_seq/logfile/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=fammi.parokkaran@students.unibe.ch

## define variables
WORKDIR="/data/users/fparokkaran/RNA_seq"
REFDIR="$WORKDIR/Reference_Genome_Data" #(GRCm39 version 110)

# Run HISAT2 to build the index for the reference genome
apptainer exec --bind /data/ \
  /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
  hisat2-build -p 4 \
  $REFDIR/Mus_musculus.GRCm39.dna.primary_assembly.fa \
  $REFDIR/Mus_musculus.GRCm39.dna.primary_assembly.index