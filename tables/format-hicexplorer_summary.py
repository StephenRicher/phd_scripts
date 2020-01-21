#!/usr/bin/env python3

import sys
import os
import pandas as pd

infile=sys.argv[1]
df = pd.read_csv(infile, sep='\t', index_col=False)

# Sort capture regions by region length and remove duplicates.
regions = (df.sort_values(by='region_length', ascending=False)
                .capture_region
                .drop_duplicates())

# Convert data to wide format for latex
df = df.pivot(index='capture_region', columns='sample',
              values='hic_pairs_per_kb')

# Reindex to order data by capture region length
df = df.reindex(list(regions))

outfile = os.path.splitext(infile)[0] + '-latex.csv'
df.to_csv(outfile)
