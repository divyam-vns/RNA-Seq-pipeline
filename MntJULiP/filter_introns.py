import pandas as pd

# Read the differential splicing output
df = pd.read_csv('multiway_out/intron.diff', sep='\t')

# Apply filtering criteria
filtered = df[(df['deltaPSI'].abs() > 0.1) & (df['qval'] < 0.05)]

# Save results
filtered.to_csv('filtered_introns.csv', index=False)
