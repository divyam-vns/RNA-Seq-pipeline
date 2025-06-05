# rMATS Event-Level Splicing Analysis Pipeline

This README describes step-by-step instructions to perform differential splicing analysis at the event level using rMATS.

## Step 1: Prepare the Environment

1. Create a working directory for rMATS:
    ```bash
    mkdir rMATS
    cd rMATS
    ```

2. Activate your rMATS Conda environment:
    ```bash
    conda activate rmats_env
    # You may need to `module load` something first depending on your HPC
    ```

## Step 2: Prepare Required Inputs

### A. Annotation File
Download or specify a GTF annotation file. Example:
```
~/genomes/gencode.v39.annotation.gtf
```

### B. BAM File Lists
Prepare two text files:

- `control_bams.txt`: comma-separated full paths to control BAMs
- `test_bams.txt`: comma-separated full paths to test BAMs

Example (`control_bams.txt`):
```
/full/path/sample1.bam,/full/path/sample2.bam,/full/path/sample3.bam
```

## Step 3: Validate Read Length
Ensure all your BAM files were generated from reads of the same length:
```bash
zcat sample.fastq.gz | head -2
```

## Step 4: Build rMATS Command

```bash
python rmats.py   --b1 control_bams.txt   --b2 test_bams.txt   --gtf ~/genomes/gencode.v39.annotation.gtf   --od control_vs_test_output   --tmp tmp_control_vs_test   --readLength 100   --nthread 10   --libType SE   --novelSS 1   --allow-clipping   > rmats.log 2>&1
```

### Explanation of Options:
- `--b1`: BAM list for control
- `--b2`: BAM list for test
- `--gtf`: GTF annotation file
- `--od`: Output directory
- `--tmp`: Temporary files directory
- `--readLength`: Read length (e.g., 100)
- `--nthread`: Number of threads
- `--libType`: `SE` for single-end or `PE` for paired-end
- `--novelSS`: Enable novel splice site detection
- `--allow-clipping`: Enable soft clipping (if aligned with STAR)

## Step 5: Run the Analysis

```bash
bash run_rmats.sh
# OR run directly:
python rmats.py ... [same options as above]
```

## Step 6: Review Output

Check the `control_vs_test_output` directory. Key files include:

- `fromGTF.SE.txt`: Skipped exons
- `fromGTF.MXE.txt`: Mutually exclusive exons
- `fromGTF.RI.txt`: Retained introns
- `fromGTF.A3SS.txt`: Alternative 3' splice sites
- `fromGTF.A5SS.txt`: Alternative 5' splice sites

Differential analysis files (e.g., `SE.MATS.JC.txt`, `SE.MATS.JCEC.txt`) include PSI, p-values, and FDR.

---