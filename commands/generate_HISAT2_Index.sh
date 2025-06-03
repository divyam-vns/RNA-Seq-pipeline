#!/bin/bash

# Set paths
BASEDIR=/path/to/project
DATADIR=${BASEDIR}/data
HISAT2_BIN=/path/to/hisat2
GENOME_INDEX_PREFIX=${BASEDIR}/hisat2_index/genome

# Create output directory
mkdir -p $(dirname $GENOME_INDEX_PREFIX)

# Build HISAT2 genome index
$HISAT2_BIN/hisat2-build ${DATADIR}/genome.fa $GENOME_INDEX_PREFIX
