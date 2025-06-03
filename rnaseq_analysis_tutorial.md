
# 🧬 RNA-seq Analysis Tutorial for Beginners

## 📁 Step 1: Set Up Your Environment

1. **Clone the Tutorial Repository:**
```bash
git clone https://github.com/griffithlab/rnaseq_tutorial.git
cd rnaseq_tutorial
```

2. **Install Required Tools:**
Ensure you have the following tools installed:
- FastQC
- HISAT2
- StringTie
- SAMtools
- DESeq2 (R package)

You can install them using `conda`, `apt`, or from official websites.

## 📥 Step 2: Download Test Data

1. **Create a Data Directory:**
```bash
mkdir -p data
cd data
```

2. **Download and Extract Test Data:**
```bash
wget http://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar
tar -xvf HBR_UHR_ERCC_ds_5pc.tar
cd ..
```

## 🔍 Step 3: Perform Quality Control

1. **Run FastQC on All FASTQ Files:**
```bash
mkdir -p qc_reports
for file in data/*.fastq.gz
do
    fastqc $file -o qc_reports/
done
```

2. **Review the Reports:**
Open the HTML reports in the `qc_reports` directory.

## 🧬 Step 4: Align Reads to Reference Genome

1. **Download Reference Genome and Annotation:**
```bash
mkdir -p reference
cd reference
wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.22.fa.gz
wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
gunzip *.gz
cd ..
```

2. **Build HISAT2 Index:**
```bash
hisat2-build reference/Homo_sapiens.GRCh37.75.dna.chromosome.22.fa reference/chr22_index
```

3. **Align Reads:**
```bash
mkdir -p alignments
for sample in UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22 \
              UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22 \
              UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22 \
              HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22 \
              HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22 \
              HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22
do
    hisat2 -x reference/chr22_index \
           -1 data/${sample}.read1.fastq.gz \
           -2 data/${sample}.read2.fastq.gz \
           -S alignments/${sample}.sam
done
```

4. **Convert SAM to BAM and Sort:**
```bash
for sam_file in alignments/*.sam
do
    bam_file=${sam_file%.sam}.bam
    sorted_bam_file=${sam_file%.sam}_sorted.bam
    samtools view -bS $sam_file > $bam_file
    samtools sort $bam_file -o $sorted_bam_file
    samtools index $sorted_bam_file
done
```

## 📊 Step 5: Quantify Gene Expression

1. **Run StringTie for Transcript Assembly and Quantification:**
```bash
mkdir -p expression
for bam_file in alignments/*_sorted.bam
do
    sample=$(basename $bam_file _sorted.bam)
    stringtie $bam_file -G reference/Homo_sapiens.GRCh37.75.gtf \
              -e -B -o expression/${sample}.gtf
done
```

## 📈 Step 6: Differential Expression Analysis

1. **Prepare for DESeq2 Analysis:**
```bash
python prepDE.py -i expression/
```

2. **Run DESeq2 in R:**
```R
library(DESeq2)
countData <- read.csv("gene_count_matrix.csv", row.names=1)
colData <- data.frame(
    row.names=colnames(countData),
    condition=c("UHR", "UHR", "UHR", "HBR", "HBR", "HBR")
)
dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~condition)
dds <- DESeq(dds)
res <- results(dds)
head(res)
```

## 💾 Step 7: Save and Upload Results to GitHub

1. **Package Your Results:**
```bash
tar -czvf rnaseq_results.tar.gz qc_reports/ alignments/ expression/ gene_count_matrix.csv
```

2. **Upload to GitHub:**
```bash
git clone https://github.com/yourusername/your_repository.git
cd your_repository
mv ../rnaseq_results.tar.gz .
git add rnaseq_results.tar.gz
git commit -m "Add RNA-seq analysis results"
git push origin main
```

---

This guide is based on the Griffith Lab’s [rnaseq_tutorial](https://github.com/griffithlab/rnaseq_tutorial).
