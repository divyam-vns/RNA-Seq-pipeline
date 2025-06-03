import pandas as pd
import re

# To add gene name
# Load GTF file
gtf_file = "/path/to/annotation.gtf"
gtf = pd.read_csv(gtf_file, sep='\t', comment='#', header=None)

# Extract gene_id and gene_name
gtf['gene_id'] = gtf[8].str.extract('gene_id "([^"]+)"')
gtf['gene_name'] = gtf[8].str.extract('gene_name "([^"]+)"')

# Create mapping
gene_mapping = gtf[['gene_id', 'gene_name']].drop_duplicates().set_index('gene_id')

# Load DESeq2 results
res = pd.read_csv("deseq2_results.csv", index_col=0)

# Add gene names
res['gene_name'] = res.index.map(lambda x: gene_mapping.loc[x.split('.')[0], 'gene_name'] if x.split('.')[0] in gene_mapping.index else 'NA')

# Save annotated results
res.to_csv("deseq2_results_annotated.csv")
