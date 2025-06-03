#!/bin/bash

# PsiCLASS performs simultaneous transcript assembly across multiple samples, enhancing the accuracy of transcript identification. 
# Set paths
BASEDIR=/path/to/project
BAM_DIR=${BASEDIR}/aligned_bams
OUTPUT_DIR=${BASEDIR}/psiclass_output
ANNOTATION_GTF=${BASEDIR}/data/annotation.gtf
PSICLASS_BIN=/path/to/psiclass

# Create output directory
mkdir -p $OUTPUT_DIR

# Generate list of BAM files
find $BAM_DIR -name "*.bam" > ${OUTPUT_DIR}/bam_list.txt

# Run PsiCLASS
$PSICLASS_BIN --lb ${OUTPUT_DIR}/bam_list.txt -o ${OUTPUT_DIR}/psiclass -p 4

# Add gene names
$PSICLASS_BIN/add-genename $ANNOTATION_GTF ${OUTPUT_DIR}/psiclass_gtf.list -o ${OUTPUT_DIR}/WithGeneNames

# Add strand information
$PSICLASS_BIN/add-strand $ANNOTATION_GTF -r < ${OUTPUT_DIR}/WithGeneNames/psiclass_vote.gtf > ${OUTPUT_DIR}/WithGeneNames/psiclass_vote.withStrand.gtf
