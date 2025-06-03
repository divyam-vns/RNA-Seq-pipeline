# Visualization of DESeq2 Differential Expression Results

This project includes an R script and a shell script to visualize results from DESeq2 differential expression analysis. The visualizations include:

- **MA Plot**
- **Volcano Plot**
- **Heatmap of Top 30 DE Genes**

## Project Structure

```
project/
├── deseq2_results.csv        # DESeq2 results with columns: log2FoldChange, padj
├── counts.csv                # Raw or normalized counts matrix (rows = genes, columns = samples)
├── sample_info.csv           # Metadata for samples (e.g., condition, batch)
├── visualize_deseq2.R        # R script for visualization
├── run_visualization.sh      # Shell script to run the R script
└── plots/                    # Output PDF plots
```

## How to Use

### 1. Make the shell script executable

```bash
chmod +x run_visualization.sh
```

### 2. Run the visualization

```bash
./run_visualization.sh
```

### 3. Check Outputs

Generated plots will be saved in the `plots/` directory as:

- `MA_plot.pdf`
- `Volcano_plot.pdf`
- `Heatmap_top30.pdf`

## Dependencies

Make sure you have the following R packages installed:

- `ggplot2`
- `pheatmap`
- `DESeq2`
- `RColorBrewer`

Install in R using:

```r
install.packages("ggplot2")
install.packages("pheatmap")
install.packages("RColorBrewer")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
```

## Contact

For questions or improvements, feel free to contribute or raise an issue.
