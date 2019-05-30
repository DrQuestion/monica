import pandas as pd
from monica.plots import barplot

n_alignment_df = pd.read_csv('/data/aalbaneseData/monica_data/output_Xylella01/monica.dataframe', index_col=(0, 1))
r_alignment_df = pd.read_csv('/data/aalbaneseData/monica_data/output_Xylella01/raw_monica.dataframe', index_col=(0, 1))
nbt=barplot._by_taxunit(n_alignment_df)
print(len(nbt))
rbt=barplot._by_taxunit(r_alignment_df)
print(len(rbt))
nbt=barplot.filter_low_reads(nbt, rbt, 3)
print(len(nbt))

