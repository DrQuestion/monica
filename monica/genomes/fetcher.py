import pandas as pd
from ete3 import NCBITaxa

import monica.genomes.tables as tables

ncbi=NCBITaxa()
PARENTS=['Fungi','Oomycota','Bacteria','Archaea','Viruses','Viroids','Nematodes','Rhizaria','Alveolata','Heterokonta']

def descendants_taxid_finder(species=[]):
    descendants = []
    for specie in species:
        for i in ncbi.get_descendant_taxa(specie):
            descendants.append(str(i))
    taxids = pd.DataFrame(descendants, columns=['taxid'])
    return taxids

def ftp_selector(mode=None, species=[]):

    species_name = []

    if not mode in ['single', 'all', 'overnight']:
        raise Exception('No mode or wrong mode specified. Please select one among the following: single, all or overnight')

    elif mode=='overnight':
        print ('Activated overnight mode')
        taxids=descendants_taxid_finder(PARENTS)
        table = tables.importer(which='genbank')
        merged_table = table.merge(taxids, on='taxid')
        for name in merged_table.loc[:, 'organism_name']:
            species_name.append(name.split(sep=' ')[0])
        merged_table['genera'] = species_name
        merged_table = merged_table.drop_duplicates(subset=['genera'], keep='last')

    elif not species:
        # only when specifically wanted single or all mode, but no species inserted
        raise Exception('You did not specify any specie. Wish to run overnight mode?')

    elif mode=='all':
        print('Activated all mode')
        taxids = descendants_taxid_finder(species)
        table=tables.importer(which='genbank')
        merged_table=table.merge(taxids, on='taxid')
        for name in merged_table.loc[:, 'organism_name']:
            species_name.append(
                ' '.join([name.split(sep=' ')[0], name.split(sep=' ')[1]]))
        merged_table['species_name'] = species_name

    elif mode=='single':
        print ('Activated single mode')
        taxids=descendants_taxid_finder(species)
        table=tables.importer(which='refseq')
        merged_table = table.merge(taxids, on='taxid')
        for name in merged_table.loc[:, 'organism_name']:
            species_name.append(
                   ' '.join([name.split(sep=' ')[0], name.split(sep=' ')[1]]))
        merged_table['species_name']=species_name
        merged_table=merged_table.drop_duplicates(subset=['species_name'], keep='last')

    merged_table.loc[:, 'ftp_path']=merged_table.loc[:, 'ftp_path']+'_genomic.fna.gz'
    return merged_table





def fetcher():
    pass

def function():
    print('lol')
