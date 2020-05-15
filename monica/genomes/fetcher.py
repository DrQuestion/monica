import os
import time
import datetime as dt

import wget
import pandas as pd
from ete3 import NCBITaxa

from . import tables


PARENTS = ['Fungi', 'Oomycota', 'Bacteria', 'Archaea', 'Viruses', 'Viroids',
           'Nematodes', 'Rhizaria', 'Alveolata', 'Heterokonta']

with open(os.path.join(os.path.join(os.path.expanduser('~'), '.monica'), '.root'), 'r') as root:
    MONICA_ROOT = root.readline()

GENOMES_PATH = os.path.join(MONICA_ROOT, 'genomes')
OLDIES_PATH = os.path.join(GENOMES_PATH, 'oldies')
EXCEPTIONS_PATH = os.path.join(GENOMES_PATH, 'exceptions')

NCBI_TAXA_UPDATE_LOG = 'ncbi_taxa_update_log'
NCBI_TAXA_DAYS_THRESHOLD = 14


def descendants_taxid_finder(species=[], focus=False):
    ncbi = NCBITaxa()
    if not ncbi_taxa_updated():
        with open(os.path.join(MONICA_ROOT, NCBI_TAXA_UPDATE_LOG), 'w+') as log:
            log.write(str(dt.date.today()))
        ncbi.update_taxonomy_database()
    descendants = []
    for specie in species:
        for i in ncbi.get_name_translator([specie])[specie]:
            descendants.append(str(i))
        for i in ncbi.get_descendant_taxa(specie):
            descendants.append(str(i))
    if not focus:
        taxids = pd.DataFrame(descendants, columns=['taxid'])
    else:
        taxids = pd.DataFrame(descendants, columns=['species_taxid'])
    return taxids


def ftp_selector(mode=None, species=[]):
    species_name = []
    ftp_path_list = []

    if mode == 'overnight':
        print('Activated overnight mode')
        genera = []
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

    elif mode == 'all':
        print('Activated all mode')
        taxids = descendants_taxid_finder(species)
        table = tables.importer(which='genbank')
        merged_table = table.merge(taxids, on='taxid')
        for name in merged_table.loc[:, 'organism_name']:
            species_name.append(
                '_'.join(name.split(sep=' ')))
        merged_table['species_name'] = species_name

    elif mode == 'single':
        print('Activated single mode')
        taxids = descendants_taxid_finder(species)
        table = tables.importer(which='refseq')
        merged_table = table.merge(taxids, on='taxid')
        for name in merged_table.loc[:, 'organism_name']:
            species_name.append(
                   '_'.join([name.split(sep=' ')[0], name.split(sep=' ')[1]]))
        merged_table['species_name'] = species_name
        merged_table = merged_table.drop_duplicates(subset=['species_name'], keep='last')

    elif mode == 'focus':
        print('Activated focus mode')
        taxids = descendants_taxid_finder(species, focus=True)
        table = tables.importer(which='genbank')
        merged_table = table.merge(taxids, on='species_taxid')

        for name, strain in zip(merged_table.loc[:, 'organism_name'], merged_table.loc[:, 'infraspecific_name']):
            if not isinstance(strain, float):
                strain = strain.split(sep='=')[1]
                if not name.endswith(strain):
                    name = name.replace('.', '')+' '+strain
                else:
                    name = name.replace(strain, '')
                    name = name.replace('.', '') + strain
            name = name.replace(' ', '_')
            species_name.append(name)

        merged_table['species_name'] = species_name
        merged_table = merged_table.drop_duplicates(subset=['species_name'], keep='last')

    # modify ftps to obtain genomic dna file
    for ftp in merged_table.loc[:, 'ftp_path']:
        filename = ftp.split(sep='/')[-1]+'_genomic.fna.gz'
        ftp_path_list.append('/'.join([ftp, filename]))
    merged_table.loc[:, 'ftp_path'] = ftp_path_list

    return merged_table


def fetcher(table, oldies_path=OLDIES_PATH, keep_genomes=None, format_genomes=None):
    t0 = time.time()

    new_genomes = []
    genomes = []

    if not os.path.exists(GENOMES_PATH):
        os.mkdir(GENOMES_PATH)
    if not os.path.exists(oldies_path):
        os.mkdir(oldies_path)
    if not os.path.exists(EXCEPTIONS_PATH):
        os.mkdir(EXCEPTIONS_PATH)

    old = os.listdir(oldies_path)

    print('Started genomes retrieval')
    if not old and not format_genomes:
        for row in table.iterrows():
            ftp = row[1]['ftp_path']
            new_header_components = [row[1][-1], row[1]['# assembly_accession'].split(sep='_')[-1]]
            if keep_genomes:
                new_filename = os.path.join(oldies_path, '_'.join(new_header_components)+'.fna.gz')
            else:
                new_filename = os.path.join(GENOMES_PATH, '_'.join(new_header_components)+'.fna.gz')
            try:
                wget.download(ftp, out=new_filename)
                genomes.append((new_filename, new_header_components))
            except FileNotFoundError:
                print('{} failed download'.format(ftp))

    elif not old and format_genomes:
        genomes_to_format = [file for file in os.listdir(format_genomes) if file.endswith('fna.gz')]
        for row in table.iterrows():
            ftp = row[1]['ftp_path']
            filename = ftp.split(sep='/')[-1]
            new_header_components = [row[1][-1], row[1]['# assembly_accession'].split(sep='_')[-1]]

            if filename in genomes_to_format:
                filename = os.path.join(format_genomes, filename)
                genomes.append((filename, new_header_components))

            else:
                if keep_genomes:
                    new_filename = os.path.join(oldies_path, '_'.join(new_header_components) + '.fna.gz')
                else:
                    new_filename = os.path.join(GENOMES_PATH, '_'.join(new_header_components) + '.fna.gz')
                try:
                    wget.download(ftp, out=new_filename)
                    genomes.append((new_filename, new_header_components))
                except FileNotFoundError:
                    print('{} failed download'.format(ftp))

    elif old and not format_genomes:
        for row in table.iterrows():
            ftp=row[1]['ftp_path']
            new_header_components=[row[1][-1], row[1]['# assembly_accession'].split(sep='_')[-1]]
            base_new_filename='_'.join(new_header_components)+'.fna.gz'

            if base_new_filename in old:
                new_filename=os.path.join(oldies_path, base_new_filename)
                genomes.append((new_filename, new_header_components))
                old.remove(base_new_filename)

            else:
                if keep_genomes:
                    new_filename = os.path.join(oldies_path, base_new_filename)
                else:
                    new_filename = os.path.join(GENOMES_PATH, base_new_filename)
                try:
                    wget.download(ftp, out=new_filename)
                    new_genomes.append(base_new_filename.split(sep='.')[0])
                    genomes.append((new_filename, new_header_components))
                except FileNotFoundError:
                    print('{} failed download'.format(ftp))
        oldies_cleaner(new_genomes, old, oldies_path)

    else:
        genomes_to_format = [file for file in os.listdir(format_genomes) if file.endswith('fna.gz')]
        for row in table.iterrows():
            ftp = row[1]['ftp_path']
            filename = ftp.split(sep='/')[-1]
            new_header_components = [row[1][-1], row[1]['# assembly_accession'].split(sep='_')[-1]]
            base_new_filename = '_'.join(new_header_components) + '.fna.gz'

            if base_new_filename in old:
                new_filename=os.path.join(oldies_path, base_new_filename)
                genomes.append((new_filename, new_header_components))
                old.remove(base_new_filename)

            else:
                if filename in genomes_to_format:
                    filename = os.path.join(format_genomes, filename)
                    genomes.append((filename, new_header_components))

                else:
                    if keep_genomes:
                        new_filename = os.path.join(oldies_path, base_new_filename)
                    else:
                        new_filename = os.path.join(GENOMES_PATH, base_new_filename)
                    try:
                        wget.download(ftp, out=new_filename)
                        new_genomes.append(base_new_filename.split(sep='.')[0])
                        genomes.append((new_filename, new_header_components))
                    except FileNotFoundError:
                        print('{} failed download'.format(ftp))
        oldies_cleaner(new_genomes, old, oldies_path)

    print('Finished genomes retrieval in {} seconds'.format(time.time()-t0))
    return genomes


def focus_fetcher(table, oldies_path=OLDIES_PATH, keep_genomes=None):
    t0 = time.time()

    genomes = []
    new_genomes = []

    old = os.listdir(oldies_path)
    temp_genomes = [file for file in os.listdir(GENOMES_PATH) if file.endswith('fna.gz')]

    print('Started genomes to focus on retrieval')
    for row in table.iterrows():
        ftp = row[1]['ftp_path']
        organism_name = row[1][-1]
        accession = row[1]['# assembly_accession'].split(sep='_')[-1]
        base_new_filename = '_'.join(['_'.join(organism_name.split('_')[0:2]), accession]) + '.fna.gz'
        if base_new_filename in old:
            new_filename=os.path.join(oldies_path, base_new_filename)
            genomes.append((new_filename, [organism_name, accession]))
            old.remove(base_new_filename)
        elif base_new_filename in temp_genomes:
            new_filename = os.path.join(GENOMES_PATH, base_new_filename)
            genomes.append((new_filename, [organism_name, accession]))
            temp_genomes.remove(base_new_filename)
        else:
            if keep_genomes:
                new_filename = os.path.join(oldies_path, base_new_filename)
            else:
                new_filename = os.path.join(GENOMES_PATH, base_new_filename)
            try:
                wget.download(ftp, out=new_filename)
                new_genomes.append(base_new_filename.split(sep='.')[0])
                genomes.append((new_filename, [organism_name, accession]))
            except FileNotFoundError:
                print('{} failed download'.format(ftp))

    oldies_cleaner(new_genomes, old, oldies_path)
    print('Finished genomes to focus on retrieval in {} seconds'.format(time.time() - t0))
    return genomes


def ncbi_taxa_updated():
    if NCBI_TAXA_UPDATE_LOG not in os.listdir(MONICA_ROOT):
        return 0
    with open(os.path.join(MONICA_ROOT, NCBI_TAXA_UPDATE_LOG), 'r') as log:
        date = log.read()
    date = dt.datetime.strptime(date, '%Y-%m-%d')
    delta = dt.datetime.now() - date
    if delta.days > NCBI_TAXA_DAYS_THRESHOLD:
        return 0
    return 1


def oldies_cleaner(new_genomes, old, oldies_path):
    # When a genome with same accession prefix but different version is downloaded, the old version is deleted
    # from the current database
    old_no_version = list(map(lambda genome: genome.split(sep='.')[0], old))
    for genome, genome_no_version in zip(old, old_no_version):
        if genome_no_version in new_genomes:
            os.remove(os.path.join(oldies_path, genome))
            print(f'Removing {genome}, new version found')
