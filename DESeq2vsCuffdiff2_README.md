
# Differential Expression Analysis: DESeq2 vs Cuffdiff2

This project demonstrates how to determine sets of differentially expressed genes using two complementary methods: **DESeq2** and **Cuffdiff2**.

## ðObjective

To compare the performance and output of two tools â€” DESeq2 and Cuffdiff2 â€” for differential gene expression analysis, and to understand their different approaches to handling transcript and gene-level data.

---

##  Methods

### 1. **DESeq2**
- **Approach**: Gene-based
- **Input**: Raw read counts per gene
- **Normalization**: Median-of-ratios method across samples
- **Output**:
  - Normalized counts per gene per sample
  - Statistical test results for differential expression
- **Limitations**: Cannot compare expression across different genes
- **Best For**: Robust gene-level testing, especially when comparing conditions or treatments

### 2. **Cuffdiff2**
- **Approach**: Transcript-based
- **Input**: Aligned reads with known transcript annotations
- **Normalization**: FPKM/TPM values (accounts for transcript length and library size)
- **Output**:
  - Expression per transcript and gene
  - Supports comparison across both genes and samples
- **Best For**: Isoform-level resolution and cross-gene comparisons

---

##  Files

- `notes.md` - Summary of methods and differences
- `data/` â-  Input data (not included here; add your own FASTQ/GTF/count files)
- `scripts/` â- R/Python scripts for running DESeq2 and parsing Cuffdiff2 output

---

## Requirements

- R with DESeq2 installed
- Cuffdiff2 (from the Cufflinks suite)
- Aligned BAM files and annotation GTF file

---

## How to Run

### DESeq2 (R)
```r
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData, colData, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
```

### Cuffdiff2 (command line)
```bash
cuffdiff -o cuffdiff_output -L control,treatment reference.gtf sample1.bam sample2.bam
```

---

## Output

- `DESeq2_output.csv`
- `cuffdiff_output/` (FPKM, differential analysis results)

---

## ðSummary: Which Tool is Better for RNA-seq?

| Feature                | DESeq2                        | Cuffdiff2                              |
|------------------------|-------------------------------|----------------------------------------|
| Analysis Level         | Gene-based                    | Transcript-based                       |
| Normalization          | Median-of-ratios              | FPKM / TPM                             |
| Handles Isoforms       | No                            | âYes                                   |
| Read Ambiguity         | Not handled                   | âHandled statistically                 |
| Cross-gene Comparison  | No                            | âYes                                   |
| Ideal Use Case         | Gene-level comparison         | Isoform-level and gene comparisons     |

**Conclusion**: Use **DESeq2** for robust, statistically strong gene-level analysis. Use **Cuffdiff2** if transcript-level resolution or isoform-specific expression is important.

---

## References

- [DESeq2 Bioconductor Page](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- [Cuffdiff2 Documentation](http://cole-trapnell-lab.github.io/cufflinks/)
