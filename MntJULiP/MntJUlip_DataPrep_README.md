
# Intron-Level Data Preparation & Pairwise Comparison with MntJUlip

This README provides step-by-step instructions to perform intron-level differential splicing analysis using MntJUlip in a Linux environment.

## 1. Environment Setup

1. Create a working directory:
```bash
mkdir MntJUlip_analysis
cd MntJUlip_analysis
```

2. Activate Conda environment with MntJUlip:
```bash
conda activate mntjulip_env
```

3. Check and save usage information (optional):
```bash
mntjulip --help | tee mntjulip_usage.txt
```

## 2. Input Requirements

MntJUlip accepts:
- A list of BAM files or pre-extracted splice junction files (.splice)
- Optional GTF annotation file for gene mapping

## 3. Extract Splice Junctions (Preferred Method)

Use `junk` (included in MntJUlip) to extract splice junctions from BAM files.

Example script:
```bash
mkdir star  # Output directory for .splice files

for bam in /path/to/bams/*.bam; do
    base=$(basename "$bam" .bam)
    if [ ! -f star/${base}.splice ]; then
        junk "$bam" > star/${base}.splice
    fi
done
```

Each `.splice` file contains:
- Chromosome, intron start/end, strand
- Read counts: total, unique, multi-mapped
- Edit statistics

## 4. Run MntJUlip (Pairwise Comparison)

Example bash script `com.mntjulip`:
```bash
BASE_DIR=/path/to/MntJUlip_analysis
ALIGN_DIR=$BASE_DIR/star
ANNOTATION=$BASE_DIR/gencode.gtf
WORK_DIR=$BASE_DIR/output
SOFTWARE_DIR=/path/to/mntjulip_installation

mntjulip \
  --input $ALIGN_DIR/*.splice \
  --annotation $ANNOTATION \
  --output $WORK_DIR \
  --threads 8 \
  --min-count 5 \
  --raw-counts-only
```

Run the script:
```bash
bash com.mntjulip
```

## 5. Output Files

- `*.introns`: List of introns and metadata
- `*.counts`: Intron read counts per sample
- `*.DSR.tsv`: Differential splicing ratio analysis
- `*.DSA.tsv`: Differential splicing abundance analysis
- `*.signif.tsv`: Significantly changed introns
- Log/JSON files: Parameters and diagnostics

## 6. Tips

- Use `.splice` files for faster reruns with different parameters.
- Use annotation (GTF) only for assigning gene names (not creating events).
- Visualize results with IGV, sashimi plots, or R/Python plots.

---
For multi-group comparison or adding covariates, see Part 2.
