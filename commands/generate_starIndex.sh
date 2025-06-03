#!/bin/bash

# Set paths
BASEDIR=/path/to/project
DATADIR=${BASEDIR}/data
STAR_BIN=/path/to/STAR
GENOME_INDEX_DIR=${BASEDIR}/star_index

# Create output directory
mkdir -p $GENOME_INDEX_DIR

# Build STAR genome index
$STAR_BIN --runThreadN 4 \
          --runMode genomeGenerate \
          --genomeDir $GENOME_INDEX_DIR \
          --genomeFastaFiles ${DATADIR}/genome.fa \
          --sjdbGTFfile ${DATADIR}/annotation.gtf \
          --sjdbOverhang 100 \
          --genomeSAindexNbases 14
