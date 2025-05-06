#!/bin/bash
# List of SRA accession numbers
ACCESSIONS=("SRR3734796" "SRR3734797" "SRR3734798")

# Loop through each accession and download the data
for accession in "${ACCESSIONS[@]}"
do
  prefetch $accession
  fasterq-dump $accession --outdir /home/projects/data
done

