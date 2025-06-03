#!/bin/bash

# Set paths
BASEDIR=/path/to/project
BAM_DIR=${BASEDIR}/aligned_bams
OUTPUT_DIR=${BASEDIR}/cuffdiff_output
ANNOTATION_GTF=${BASEDIR}/data/annotation.gtf
CUFFDIFF_BIN=/path/to/cuffdiff

# Create output directory
mkdir -p $OUTPUT_DIR

# Define sample groups
GROUP1_BAMS="${BAM_DIR}/sample1.bam,${BAM_DIR}/sample2.bam"
GROUP2_BAMS="${BAM_DIR}/sample3.bam,${BAM_DIR}/sample4.bam"

# Run Cuffdiff
$CUFFDIFF_BIN -o $OUTPUT_DIR -p 4 $ANNOTATION_GTF $GROUP1_BAMS $GROUP2_BAMS
