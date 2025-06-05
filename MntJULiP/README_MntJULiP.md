
# Differential Splicing Analysis using MntJULiP

This guide provides step-by-step instructions for performing differential splicing analysis using MntJULiP, including pairwise and multiway comparisons with covariate-aware modeling (e.g., sex, age, BMI).

---

## 🔧 Requirements

- MntJULiP installed
- Reference annotation file (GTF/GFF3)
- BAM files for RNA-seq data
- Metadata file (with condition and covariates)

---

## ✅ Step 1: Generate Splice Files

```bash
mntjulip extract \
  -a annotation.gtf \
  -s samples.tsv \
  -o ./splice_files/ \
  -t 24
```

**Example `samples.tsv`:**
```
sample_name\tbam_path
S1\t/path/to/S1.bam
S2\t/path/to/S2.bam
```

---

## ✅ Step 2: Pairwise Comparison (No Covariates)

**`splice.list`:**
```
sample1.splice\tConditionA
sample2.splice\tConditionA
sample3.splice\tConditionB
sample4.splice\tConditionB
```

**Command:**
```bash
mntjulip dsr \
  -a annotation.gtf \
  -s splice.list \
  -o pairwise_out/ \
  -t 24
```

---

## ✅ Step 3: Pairwise Comparison with Covariates

**Enhanced `splice.list`:**
```
sample1.splice\tMonth1\tFemale
sample2.splice\tMonth1\tFemale
sample3.splice\tMonth4\tMale
sample4.splice\tMonth4\tMale
```

**Command:**
```bash
mntjulip dsr \
  -a annotation.gtf \
  -s splice.list \
  -o pairwise_cov_out/ \
  -t 24 \
  --covariates 2
```

---

## ✅ Step 4: Multiway Comparison (with/without covariates)

**Multi-condition `splice.list`:**
```
sample1.splice\tMonth1\tFemale
sample2.splice\tMonth2\tMale
sample3.splice\tMonth4\tFemale
```

**Command with covariates:**
```bash
mntjulip dsr \
  -a annotation.gtf \
  -s splice.list \
  -o multiway_out/ \
  -t 24 \
  --covariates 2
```

**Without covariates:**
```bash
mntjulip dsr \
  -a annotation.gtf \
  -s splice.list \
  -o multiway_no_cov_out/ \
  -t 24
```

---

## ✅ Step 5: Python Filtering Example

Filter introns with |ΔPSI| > 0.1 and q-value < 0.05:

```python
import pandas as pd

df = pd.read_csv('multiway_out/intron.diff', sep='\t')
filtered = df[(df['deltaPSI'].abs() > 0.1) & (df['qval'] < 0.05)]
filtered.to_csv('filtered_introns.csv', index=False)
```

---

## 🧪 Optional: Generate Splice Files for New Samples

```bash
mntjulip extract \
  -a annotation.gtf \
  -s new_samples.tsv \
  -o ./splice_files/multi/ \
  -t 24
```

---

## 📁 Files Included

- `mntjulip_analysis.sh`: Bash script to run pairwise and multiway comparisons
- `filter_introns.py`: Python script to filter significant differential introns
- `README.md`: This documentation file
