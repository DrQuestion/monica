import os
from collections import Counter

import pandas as pd
import mappy

from monica.genomes.database import DATABASE

# db='/home/drq/Desktop/genome.fna'
# query='/home/drq/Desktop/query.fa'

def indexer (database=DATABASE):
    index = mappy.Aligner(fn_idx_in=database, preset='map-ont')
    return index

def aligner (query_folder, index):
    allignment=dict()

    os.chdir(query_folder)
    samples=os.listdir('.')

    for sample in samples:
        allignment[sample]=dict()
        for name, seq, qual in mappy.fastx_read(sample):
            for hit in index.map(seq):
                tax_unit = hit.ctg.split(sep=':')[0]
                accession = hit.ctg.split(sep=':')[1]
                if tax_unit in allignment[sample]:
                    allignment[sample][tax_unit].update({accession: 1})
                else:
                    allignment[sample][tax_unit] = Counter({accession: 1})

    return(allignment)

def allignment_to_data_frame (allignment):
    data_frame = pd.concat({k: pd.DataFrame(v).unstack() for k, v in allignment.items()}, axis=1)
    return data_frame

#Implement to test if a count system based on sequenced bases works better/well too

if __name__=='__main__':
    db='/home/drq/Desktop/temp/genome.fna'
    query='/home/drq/Desktop/temp/query'
    al=aligner(query, db)
    print(al)