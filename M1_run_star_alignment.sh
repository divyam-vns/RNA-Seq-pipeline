#!/bin/bash

#===============================
# STAR Alignment Automation Script
# Includes logging and error handling
#===============================

# Set number of threads
THREADS=4

# Define paths
GENOME_DIR="/home/projects/Data/M1_exam/data_m1_exam/star"
FASTQ_DIR="/home/projects/Data/M1_exam/data_m1_exam/transcriptomics_proj1_data"
OUT_DIR="/home/projects/Data/M1_exam/Results"
LOG_FILE="$OUT_DIR/star_alignment.log"

# Ensure output directory exists
mkdir -p "$OUT_DIR"

# Initialize log file
echo "STAR Alignment Log - $(date)" > "$LOG_FILE"
echo "--------------------------------------" >> "$LOG_FILE"

# Loop through all FASTQ files
for fq in "$FASTQ_DIR"/*.fastq.gz; do
    sample=$(basename "$fq" .fastq.gz)
    timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[$timestamp] Starting alignment for $sample..." | tee -a "$LOG_FILE"

    STAR \
        --runThreadN "$THREADS" \
        --genomeDir "$GENOME_DIR" \
        --readFilesIn "$fq" \
        --readFilesCommand zcat \
        --outFileNamePrefix "$OUT_DIR/${sample}_" \
        --outSAMtype BAM SortedByCoordinate >> "$LOG_FILE" 2>&1

    if [ $? -eq 0 ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Alignment completed for $sample." | tee -a "$LOG_FILE"
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Alignment failed for $sample." | tee -a "$LOG_FILE"
    fi

    echo "--------------------------------------" >> "$LOG_FILE"
done

echo "All alignments finished. Log saved to $LOG_FILE."

