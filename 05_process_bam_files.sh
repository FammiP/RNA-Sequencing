#!/bin/bash

#================================================================
# SCRIPT NAME:
#   05_process_bam_files.sh
#
# DESCRIPTION:
#   This script converts SAM files generated from RNA-seq mapping 
#   into compressed and sorted BAM files using Samtools. It also 
#   indexes the sorted BAM files for efficient access during downstream 
#   analysis. Optionally, the intermediate unsorted BAM files are 
#   deleted to save disk space.
#
# USAGE:
#   sbatch 05_process_bam_files.sh
#
# REQUIREMENTS:
#   - Samtools available via the specified Apptainer container
#   - Access to a compute cluster with SLURM scheduler
#   - Existing SAM files generated from the mapping step
#
# INPUT:
#   - SAM files located in: $SAMDIR
#
# OUTPUT:
#   - Sorted BAM files in: $BAMDIR
#   - BAM index files (.bai) in: $BAMDIR
#
# NOTES:
#   - Ensure Samtools is available via the container path.
#   - Update paths and email address if necessary.
#================================================================


#SBATCH --job-name=SAM_to_BAM
#SBATCH --time=120:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/fparokkaran/RNA_seq/logfile/%x_%j.out
#SBATCH --error=/data/users/fparokkaran/RNA_seq/logfile/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=fammi.parokkaran@students.unibe.ch

# Define paths
WORKDIR="/data/users/fparokkaran/RNA_seq"
SAMDIR="$WORKDIR/02_Mapping/SAM_files"
BAMDIR="$WORKDIR/02_Mapping/BAM_Files"
CONTAINER="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"
THREADS=8
MEM="4G"
# Create output directory if not already existing
mkdir -p "$BAMDIR"

for SAMFILE in "$SAMDIR"/*_mappedReads.sam; do
    BASENAME=$(basename "$SAMFILE" _mappedReads.sam)

    # Convert SAM to BAM
    BAMFILE="$BAMDIR/${BASENAME}_mappedReads.bam"
    apptainer exec --bind /data/ "$CONTAINER" samtools view -hbS "$SAMFILE" > "$BAMFILE"

    # Sort BAM File
    SORTED_BAMFILE="$BAMDIR/${BASENAME}_sorted.bam"
    apptainer exec --bind /data/ "$CONTAINER" samtools sort -m $MEM -@ $THREADS \
        -o "$SORTED_BAMFILE" "$BAMFILE"

    # Index Sorted BAM File
    echo "Indexing $SORTED_BAMFILE..."
    apptainer exec --bind /data/ "$CONTAINER" samtools index "$SORTED_BAMFILE"

    # Optional: Clean up unsorted BAM files to save space
    echo "Cleaning up $BAMFILE..."
    rm "$BAMFILE"
done

