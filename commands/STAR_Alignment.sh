#!/bin/bash

# Set paths
BASEDIR=/path/to/project
DATADIR=${BASEDIR}/data
STAR_BIN=/path/to/STAR
GENOME_INDEX_DIR=${BASEDIR}/star_index
OUTDIR=${BASEDIR}/star_output

# Create output directory
mkdir -p $OUTDIR

# Sample list
SAMPLES=("sample1" "sample2" "sample3")

# Run STAR alignment
for SAMPLE in "${SAMPLES[@]}"; do
  $STAR_BIN --runThreadN 4 \
            --genomeDir $GENOME_INDEX_DIR \
            --readFilesIn ${DATADIR}/${SAMPLE}.fastq.gz \
            --readFilesCommand zcat \
            --outFileNamePrefix ${OUTDIR}/${SAMPLE}_ \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMstrandField intronMotif
done
