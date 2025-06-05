#!/bin/bash

############################################
# MntJULiP Differential Splicing Pipeline  #
# Supports:                                #
# - Pairwise comparison                    #
# - Multiway comparison                    #
# - Covariate-aware modeling (sex, age)    #
############################################

# SET VARIABLES
ANNOTATION="annotation.gtf"
THREADS=24

# Step 1: Extract splice files
# samples.tsv: 2-column file with sample name and BAM path
echo "[Step 1] Extracting splice files..."
mntjulip extract \
  -a $ANNOTATION \
  -s samples.tsv \
  -o splice_files/ \
  -t $THREADS

# Step 2: Pairwise comparison without covariates
# splice_pairwise.list: 2 columns: splice path and condition (e.g., Month1 or Month4)
echo "[Step 2] Running pairwise comparison (no covariates)..."
mntjulip dsr \
  -a $ANNOTATION \
  -s splice_pairwise.list \
  -o results/pairwise_no_cov/ \
  -t $THREADS

# Step 3: Pairwise comparison with covariates (e.g., sex)
# splice_pairwise_cov.list: 3 columns: splice path, condition, covariate (e.g., sex)
echo "[Step 3] Running pairwise comparison (with covariates)..."
mntjulip dsr \
  -a $ANNOTATION \
  -s splice_pairwise_cov.list \
  -o results/pairwise_with_cov/ \
  -t $THREADS \
  --covariates 2

# Step 4: Multiway comparison without covariates
# splice_multi.list: 2 columns: splice path and condition (e.g., Month1, Month2, Month4)
echo "[Step 4] Running multiway comparison (no covariates)..."
mntjulip dsr \
  -a $ANNOTATION \
  -s splice_multi.list \
  -o results/multiway_no_cov/ \
  -t $THREADS

# Step 5: Multiway comparison with covariates (e.g., sex, age, BMI)
# splice_multi_cov.list: >2 columns: splice path, condition, covariate1, covariate2, etc.
echo "[Step 5] Running multiway comparison (with covariates)..."
mntjulip dsr \
  -a $ANNOTATION \
  -s splice_multi_cov.list \
  -o results/multiway_with_cov/ \
  -t $THREADS \
  --covariates 2

# Step 6: Optional filtering of results (example Python snippet)
echo "[INFO] You can filter results using the following Python script:"
echo ""
echo "import pandas as pd"
echo "df = pd.read_csv('results/multiway_with_cov/intron.diff', sep='\t')"
echo "filtered = df[(df['deltaPSI'].abs() > 0.1) & (df['qval'] < 0.05)]"
echo "filtered.to_csv('filtered_introns.csv', index=False)"
echo ""
echo "[Done]"
