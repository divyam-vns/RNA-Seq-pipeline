
# RNA-seq Differential Expression Analysis Pipeline

## Overview
This pipeline provides a clear, customizable workflow for RNA-seq differential expression analysis using:
- STAR and HISAT2 for read alignment
- PsiCLASS for transcript assembly
- Cuffdiff and DESeq2 for differential expression analysis

---

## 1. Quality Control
(Optional but recommended)

```bash
# Run FastQC
fastqc /path/to/raw_reads/*.fastq.gz -o /path/to/qc_reports
```

---

## 2. Reference Genome Indexing

### STAR Indexing

```bash
STAR --runThreadN 4      --runMode genomeGenerate      --genomeDir /path/to/star_index      --genomeFastaFiles /path/to/genome.fa      --sjdbGTFfile /path/to/annotation.gtf      --sjdbOverhang 100
```

### HISAT2 Indexing

```bash
hisat2-build /path/to/genome.fa /path/to/hisat2_index/genome
```

---

## 3. Read Alignment

### STAR

```bash
STAR --runThreadN 4      --genomeDir /path/to/star_index      --readFilesIn sample.fastq.gz      --readFilesCommand zcat      --outSAMtype BAM SortedByCoordinate
```

### HISAT2

```bash
hisat2 -p 4 -x /path/to/hisat2_index/genome        -U sample.fastq.gz        -S output.sam
```

---

## 4. Transcript Assembly (PsiCLASS)

```bash
psiclass --lb bam_list.txt -o output_dir/psiclass -p 4
psiclass/add-genename annotation.gtf psiclass_gtf.list -o output_dir/WithGeneNames
psiclass/add-strand annotation.gtf -r < input.gtf > output.gtf
```

---

## 5. Differential Expression

### Cuffdiff

```bash
cuffdiff -o cuffdiff_output -p 4 annotation.gtf group1.bam group2.bam
```

### DESeq2 (in R)

```r
library(DESeq2)
dds <- DESeqDataSet(se, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
write.csv(as.data.frame(res), "deseq2_results.csv")
```

---

## 6. Annotate Results

```python
# Annotate DESeq2 results with gene names
import pandas as pd

gtf = pd.read_csv("annotation.gtf", sep='\t', comment='#', header=None)
gtf['gene_id'] = gtf[8].str.extract('gene_id "([^"]+)"')
gtf['gene_name'] = gtf[8].str.extract('gene_name "([^"]+)"')
gene_mapping = gtf[['gene_id', 'gene_name']].drop_duplicates().set_index('gene_id')

res = pd.read_csv("deseq2_results.csv", index_col=0)
res['gene_name'] = res.index.map(lambda x: gene_mapping.loc[x.split('.')[0], 'gene_name'] if x.split('.')[0] in gene_mapping.index else 'NA')
res.to_csv("deseq2_results_annotated.csv")
```

---

## Notes
- Replace paths and filenames with actual values relevant to your dataset.
- Ensure required tools (STAR, HISAT2, PsiCLASS, Cuffdiff, R with DESeq2) are installed and accessible in your PATH.
