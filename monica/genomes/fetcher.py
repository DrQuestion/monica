import pandas as pd
from ete3 import NCBITaxa

import monica.genomes.tables as tables

ncbi=NCBITaxa()
PARENTS=['Fungi','Oomycota','Bacteria','Archaea','Viruses','Viroids','Nematodes','Rhizaria','Alveolata','Heterokonta']

def descendants_taxid_finder(species=[]):
    descendants = []
    for specie in species:
        for i in ncbi.get_descendant_taxa(specie):
            descendants.append(i)
    taxids = pd.DataFrame(descendants, columns=['taxid'])
    return taxids


def ftp_selector(mode=None, species=[]):

    species_name = []

    if not mode in ['single', 'all', 'overnight']:
        raise Exception('No mode or wrong mode specified. Please select one among the following: single, all or overnight')

    if mode=='overnight':
        print ('Activated overnight mode')
        taxids=descendants_taxid_finder(PARENTS)
        table = tables.importer(which='genbank')
        merged_table = table.merge(taxids, on='taxid')

    elif not species:
        # only when specifically wanted single or all mode, but no species inserted
        raise Exception('You did not specify any specie. Wish to run overnight mode?')

    elif mode=='all':
        print('Activated all mode')
        taxids = descendants_taxid_finder(species)
        table=tables.importer(which='genbank')
        merged_table=table.merge(taxids, on='taxid')

    elif mode=='single':
        print ('Activated single mode')
        taxids=descendants_taxid_finder(species)
        table=tables.importer(which='refseq')
        merged_table = table.merge(taxids, on='taxid')
        for row in merged_table.iterrows():
            species_name.append(
                ' '.join([row['organism_name'].split(sep=' ')[0], row['organism_name'].split(sep=' ')[1]]))




def fetcher():
    pass

def function():
    print('lol')
