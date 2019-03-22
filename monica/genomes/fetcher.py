import pandas as pd
from ete3 import NCBITaxa

import tables.importer

ncbi=NCBITaxa()

def ftp_selector(mode=None, species=[]):
    descendants = []
    species_name = []
    for specie in species:
        for i in ncbi.get_descendant_taxa(specie):
            descendants.append(i)
    taxids = pd.DataFrame(descendants, columns=['taxid'])

    #return descendants

    if not mode:
        table=tables.importer(which='refseq')

        merged_table = table.merge(taxids, on='taxid')

        for row in merged_table.iterrows():
            species_name.append(
                ' '.join([row['organism_name'].split(sep=' ')[0], row['organism_name'].split(sep=' ')[1]]))
    else:
        table=tables.importer(which='genbank')

    merged_table=table.merge(taxids, on='taxid')

    # if mode=='all':
    #     pass
    # elif mode=='overnight':
    #     pass
    # else:
    #     pass



def fetcher():
    pass

def function():
    print('lol')
