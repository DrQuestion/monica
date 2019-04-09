import os
import gzip
import time
import datetime as dt

import wget
import pandas as pd
from ete3 import NCBITaxa
from Bio import SeqIO

import monica.genomes.tables as tables


PARENTS=['Fungi','Oomycota','Bacteria','Archaea','Viruses','Viroids','Nematodes','Rhizaria','Alveolata','Heterokonta']
GENOMES_PATH=os.path.join(os.path.dirname(__file__), 'genomes')
OLDIES_PATH=os.path.join(GENOMES_PATH, 'oldies')
CWD=os.getcwd()

NCBI_TAXA_UPDATE_LOG='ncbi_taxa_update_log'
NCBI_TAXA_DAYS_THRESHOLD=14


def descendants_taxid_finder(species=[]):
    ncbi = NCBITaxa()
    if not ncbi_taxa_updated():
        with open(os.path.join(os.path.dirname(__file__), NCBI_TAXA_UPDATE_LOG), 'w+') as log:
            log.write(str(dt.date.today()))
        ncbi.update_taxonomy_database()
    descendants = []
    for specie in species:
        for i in ncbi.get_descendant_taxa(specie):
            descendants.append(str(i))
    taxids = pd.DataFrame(descendants, columns=['taxid'])
    return taxids


def ftp_selector(mode=None, species=[]):
    species_name = []
    ftp_path_list = []

    if not mode in ['single', 'all', 'overnight']:
        raise Exception('No mode or wrong mode specified. Please select one among the following: single, all or overnight')

    elif mode=='overnight':
        print ('Activated overnight mode')
        genera=[]
        taxids=descendants_taxid_finder(PARENTS)
        table = tables.importer(which='refseq')
        merged_table = table.merge(taxids, on='taxid')
        for name in merged_table.loc[:, 'organism_name']:
            splitted_name=name.split(sep=' ')
            genera.append(splitted_name[0])
            species_name.append('_'.join([splitted_name[0],splitted_name[1]]))
        merged_table['genera'] = genera
        merged_table['species_name']=species_name
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
                '_'.join([name.split(sep=' ')[0], name.split(sep=' ')[1]]))
        merged_table['species_name'] = species_name

    elif mode=='single':
        print ('Activated single mode')
        taxids=descendants_taxid_finder(species)
        table=tables.importer(which='refseq')
        merged_table = table.merge(taxids, on='taxid')
        for name in merged_table.loc[:, 'organism_name']:
            species_name.append(
                   '_'.join([name.split(sep=' ')[0], name.split(sep=' ')[1]]))
        merged_table['species_name']=species_name
        merged_table=merged_table.drop_duplicates(subset=['species_name'], keep='last')
        # TODO problem with keeping last? Xylella fastidiosa genome (003352785.1) too short

    # modify ftps to obtain genomic dna file
    for ftp in merged_table.loc[:, 'ftp_path']:
        filename=ftp.split(sep='/')[-1]+'_genomic.fna.gz'
        ftp_path_list.append('/'.join([ftp,filename]))
    merged_table.loc[:, 'ftp_path']=ftp_path_list

    return merged_table


def fetcher(table):
    oldies=[]
    new_genomes=[]

    if not os.path.exists(GENOMES_PATH):
        os.mkdir(GENOMES_PATH)
        os.mkdir(OLDIES_PATH)

    os.chdir(GENOMES_PATH)
    old=os.listdir(OLDIES_PATH)

    if not old:
        for row in table.iterrows():
            ftp=row[1]['ftp_path']
            filename=ftp.split(sep='/')[-1]
            new_header_components=[row[1][-1], row[1]['# assembly_accession'].split(sep='_')[-1]]
            new_filename='_'.join(new_header_components)+'.fna.gz'
            new_header = ':'.join(new_header_components)
            print(f'Started {filename} download')
            t0 = time.time()
            wget.download(ftp)
            header_modifier(filename, new_filename, new_header)

    else:
        for row in table.iterrows():
            ftp=row[1]['ftp_path']
            filename=ftp.split(sep='/')[-1]
            new_header_components=[row[1][-1], row[1]['# assembly_accession'].split(sep='_')[-1]]
            new_filename='_'.join(new_header_components)+'.fna.gz'

            if new_filename in old:
                oldies.append(new_filename)
                old.remove(new_filename)
                print(f'{filename} already present in oldies as {new_filename}')

            else:
                new_header = ':'.join(new_header_components)
                print(f'Started {filename} download')
                wget.download(ftp)
                header_modifier(filename, new_filename, new_header)
                new_genomes.append(new_filename.split(sep='.')[0])

    print(f'Finished genomes retrieval')
    os.chdir(CWD)
    oldies_cleaner(new_genomes, old)
    return oldies


def header_modifier(file, new_filename, new_header):
    with gzip.open(file, 'rt') as old, gzip.open(new_filename, 'wt') as new:
        for seq_record in SeqIO.parse(old, 'fasta'):
            seq_record.id=new_header
            SeqIO.write(seq_record, new, 'fasta')
    os.remove(file)


def ncbi_taxa_updated():
    if not NCBI_TAXA_UPDATE_LOG in os.listdir(os.path.dirname(__file__)):
        return 0
    with open (os.path.join(os.path.dirname(__file__), NCBI_TAXA_UPDATE_LOG), 'r') as log:
        date = log.read()
    date = dt.datetime.strptime(date, '%Y-%m-%d')
    delta = dt.datetime.now() - date
    if delta.days > NCBI_TAXA_DAYS_THRESHOLD:
        return 0
    return 1


def oldies_cleaner(new_genomes, old):
    # When a genome with same accession prefix but different version is downloaded, the old version is deleted
    # from the current database
    old_no_version=list(map(lambda genome: genome.split(sep='.')[0], old))
    for genome, genome_no_version in zip(old, old_no_version):
        if genome_no_version in new_genomes:
            os.remove(os.path.join(OLDIES_PATH, genome))
            print(f'Removing {genome}, new version found')
