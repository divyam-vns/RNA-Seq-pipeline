# RNA-Seq STAR Alignment Pipeline

This repository contains a shell script to automate the alignment of RNA-seq reads using STAR. It includes logging, error handling, and organizes output BAM files for downstream analysis.

## 📁 Directory Structure

Before running the script, ensure your directory is organized as follows:

```
/home/projects/Data/M1_exam/
├── data_m1_exam/
│   ├── transcriptomics_proj1_data/
│   │   ├── 1mo_Rep1.fastq.gz
│   │   ├── ...
│   │   ├── region.fa
│   │   ├── gencode.vM30.region.gtf
│   ├── star/  # STAR genome index (generated with --runMode genomeGenerate)
├── Results/   # Output directory for aligned BAMs and logs
├── run_star_alignment.sh
```

## 🛠 Requirements

- STAR (tested with version 2.7.11b)
- Bash
- zcat
- Sufficient memory and disk space

## 🚀 Usage

1. Give the script execution permissions:

   ```bash
   chmod +x run_star_alignment.sh
   ```

2. Run the script:

   ```bash
   ./run_star_alignment.sh
   ```

3. Aligned and sorted BAM files will be saved in the Results/ directory with logs written to star_alignment.log.

## 🧪 Example Output

For each FASTQ file:
- A sorted BAM file:  Results/sampleName_Aligned.sortedByCoord.out.bam
- STAR logs:          Results/star_alignment.log

## 📓 Notes

- STAR genome index must be pre-generated in the star/ directory using:

  ```bash
  STAR --runThreadN 4 \
       --runMode genomeGenerate \
       --genomeDir star \
       --genomeFastaFiles region.fa \
       --sjdbGTFfile gencode.vM30.region.gtf \
       --sjdbOverhang 99
  ```

- If any alignment fails, the error will be captured in the log.

## 📬 Contact

For issues or suggestions, please open an issue or contact the project maintainer.
