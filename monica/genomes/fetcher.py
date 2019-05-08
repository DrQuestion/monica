import os
import gzip
# import shutil
import pickle
import time
import datetime as dt

import wget
import pandas as pd
from ete3 import NCBITaxa
from Bio import SeqIO

from . import tables


PARENTS = ['Fungi', 'Oomycota', 'Bacteria', 'Archaea', 'Viruses', 'Viroids',
           'Nematodes', 'Rhizaria', 'Alveolata', 'Heterokonta']
GENOMES_PATH = os.path.join(os.path.dirname(__file__), 'genomes')
OLDIES_PATH = os.path.join(GENOMES_PATH, 'oldies')
# EXCEPTIONS_PATH=os.path.join(GENOMES_PATH, 'exeptions')
# CWD=os.getcwd()

NCBI_TAXA_UPDATE_LOG = 'ncbi_taxa_update_log'
NCBI_TAXA_DAYS_THRESHOLD = 14


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

    if mode == 'overnight':
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

    elif mode == 'all':
        print('Activated all mode')
        taxids = descendants_taxid_finder(species)
        table=tables.importer(which='genbank')
        merged_table=table.merge(taxids, on='taxid')
        for name in merged_table.loc[:, 'organism_name']:
            species_name.append(
                '_'.join([name.split(sep=' ')[0], name.split(sep=' ')[1]]))
        merged_table['species_name'] = species_name

    elif mode == 'single':
        print ('Activated single mode')
        taxids=descendants_taxid_finder(species)
        table=tables.importer(which='refseq')
        merged_table = table.merge(taxids, on='taxid')
        for name in merged_table.loc[:, 'organism_name']:
            species_name.append(
                   '_'.join([name.split(sep=' ')[0], name.split(sep=' ')[1]]))
        merged_table['species_name']=species_name
        merged_table=merged_table.drop_duplicates(subset=['species_name'], keep='last')

    # modify ftps to obtain genomic dna file
    for ftp in merged_table.loc[:, 'ftp_path']:
        filename = ftp.split(sep='/')[-1]+'_genomic.fna.gz'
        ftp_path_list.append('/'.join([ftp, filename]))
    merged_table.loc[:, 'ftp_path'] = ftp_path_list

    return merged_table


def fetcher(table, oldies_path=OLDIES_PATH, keep_genomes=None, format_genomes=None):
    # oldies = []
    t0=time.time()

    new_genomes = []
    genomes = []

    if not os.path.exists(GENOMES_PATH):
        os.mkdir(GENOMES_PATH)
    if not os.path.exists(oldies_path):
        os.mkdir(oldies_path)
    # if not os.path.exists(EXCEPTIONS_PATH):
    #     os.mkdir(EXCEPTIONS_PATH)

    # os.chdir(GENOMES_PATH)
    old = os.listdir(oldies_path)

    print('Started genomes retrieval')
    if not old and not format_genomes:
        # current_genomes_length = dict()
        for row in table.iterrows():
            ftp=row[1]['ftp_path']
        #    filename=ftp.split(sep='/')[-1]
            new_header_components=[row[1][-1], row[1]['# assembly_accession'].split(sep='_')[-1]]
            if keep_genomes:
                new_filename = os.path.join(oldies_path, '_'.join(new_header_components)+'.fna.gz')
            else:
                new_filename = os.path.join(GENOMES_PATH, '_'.join(new_header_components)+'.fna.gz')

        #    new_header = ':'.join(new_header_components)
            wget.download(ftp, out=new_filename)
            genomes.append((new_filename, new_header_components))
        #     try:
        #         genome_length = header_modifier(filename, new_filename, new_header)
        #         current_genomes_length[new_header_components[1]] = genome_length
        #     except:
        #         print(f'{filename} (ftp:{ftp}) failed')
        #         shutil.move(filename, EXCEPTIONS_PATH)
        #
        # pickle.dump(current_genomes_length, open(os.path.join(GENOMES_PATH, 'current_genomes_length.pkl'), 'wb'))

    elif not old and format_genomes:
        # print(f'No old genomes folder provided, formatting genomes in {format_genomes} and downloading new ones')
        # current_genomes_length = dict()
        genomes_to_format = [file for file in os.listdir(format_genomes) if file.endswith('fna.gz')]
        for row in table.iterrows():
            ftp = row[1]['ftp_path']
            filename = ftp.split(sep='/')[-1]
            new_header_components = [row[1][-1], row[1]['# assembly_accession'].split(sep='_')[-1]]
        #    new_header = ':'.join(new_header_components)

            if filename in genomes_to_format:
                filename = os.path.join(format_genomes, filename)
                genomes.append((filename, new_header_components))
        #        try:
        #            filename = os.path.join(format_genomes, filename)
        #            genome_length = header_modifier(filename, new_filename, new_header, format=True)
        #            current_genomes_length[new_header_components[1]] = genome_length
        #        except:
        #            print(f'{filename} failed formatting')
        #            continue

            else:
                # print(f'Started {filename} download')
                if keep_genomes:
                    new_filename = os.path.join(oldies_path, '_'.join(new_header_components) + '.fna.gz')
                else:
                    new_filename = os.path.join(GENOMES_PATH, '_'.join(new_header_components) + '.fna.gz')
                wget.download(ftp, out=new_filename)
                genomes.append((new_filename, new_header_components))
        #        try:
        #            genome_length = header_modifier(filename, new_filename, new_header)
        #            current_genomes_length[new_header_components[1]] = genome_length
        #        except:
        #            print(f'{filename} (ftp:{ftp}) failed')
        #            shutil.move(filename, EXCEPTIONS_PATH)

        # pickle.dump(current_genomes_length, open(os.path.join(GENOMES_PATH, 'current_genomes_length.pkl'), 'wb'))

    elif old and not format_genomes:
        # print(f'Provided an old genomes folder, downloading necessary genomes not there')
        # genomes_length = pickle.load(open(os.path.join(oldies_path, 'genomes_length.pkl'), 'rb'))
        # current_genomes_length = dict()
        for row in table.iterrows():
            ftp=row[1]['ftp_path']
        #    filename=ftp.split(sep='/')[-1]
            new_header_components=[row[1][-1], row[1]['# assembly_accession'].split(sep='_')[-1]]
            base_new_filename='_'.join(new_header_components)+'.fna.gz'

            if base_new_filename in old:
                new_filename=os.path.join(oldies_path, base_new_filename)
                genomes.append((new_filename, new_header_components))
                old.remove(base_new_filename)
        #        current_genomes_length[new_header_components[1]] = genomes_length[new_header_components[1]]
                # print(f'{filename} already present in oldies as {new_filename}')

            else:
        #        new_header = ':'.join(new_header_components)
                # print(f'Started {filename} download')
                if keep_genomes:
                    new_filename = os.path.join(oldies_path, base_new_filename)
                else:
                    new_filename = os.path.join(GENOMES_PATH, base_new_filename)
                wget.download(ftp, out=new_filename)
                new_genomes.append(base_new_filename.split(sep='.')[0])
                genomes.append((new_filename, new_header_components))
        #         try:
        #             genome_length = header_modifier(filename, new_filename, new_header)
        #             new_genomes.append(base_new_filename.split(sep='.')[0])
        #             current_genomes_length[new_header_components[1]] = genome_length
        #         except:
        #             print(f'{filename} (ftp:{ftp}) failed')
        #             shutil.move(filename, EXCEPTIONS_PATH)
        #
        # pickle.dump(current_genomes_length, open(os.path.join(GENOMES_PATH, 'current_genomes_length.pkl'), 'wb'))
        oldies_cleaner(new_genomes, old, oldies_path)

    else:
        # print(f'Provided an old genomes folder,
        # if necessary formatting genomes in {format_genomes} or downloading them')
        # genomes_length = pickle.load(open(os.path.join(oldies_path, 'genomes_length.pkl'), 'rb'))
        # current_genomes_length = dict()
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
                # current_genomes_length[new_header_components[1]] = genomes_length[new_header_components[1]]
                # print(f'{filename} already present in oldies as {new_filename}')

            else:
        #        new_header = ':'.join(new_header_components)
                if filename in genomes_to_format:
                    # print(f'Formatting {filename} found in {format_genomes}')
                    filename = os.path.join(format_genomes, filename)
                    genomes.append((filename, new_header_components))
        #            try:
        #                genome_length= header_modifier(filename, new_filename, new_header, format=True)
        #                new_genomes.append(new_filename.split(sep='.')[0])
        #                current_genomes_length[new_header_components[1]] = genome_length
        #            except:
        #                print(f'{filename} failed formatting')
        #                continue

                else:
                    # print(f'Started {filename} download')
                    if keep_genomes:
                        new_filename = os.path.join(oldies_path, base_new_filename)
                    else:
                        new_filename = os.path.join(GENOMES_PATH, base_new_filename)
                    wget.download(ftp, out=new_filename)
                    new_genomes.append(base_new_filename.split(sep='.')[0])
                    genomes.append((new_filename, new_header_components))
        #             try:
        #                 genome_length = header_modifier(filename, new_filename, new_header)
        #                 new_genomes.append(new_filename.split(sep='.')[0])
        #                 current_genomes_length[new_header_components[1]] = genome_length
        #             except:
        #                 print(f'{filename} (ftp:{ftp}) failed')
        #                 shutil.move(filename, EXCEPTIONS_PATH)
        #
        # pickle.dump(current_genomes_length, open(os.path.join(GENOMES_PATH, 'current_genomes_length.pkl'), 'wb'))
        oldies_cleaner(new_genomes, old, oldies_path)

    print('Finished genomes retrieval in {} seconds'.format(time.time()-t0))
    # os.chdir(CWD)
    return genomes


def header_modifier(filename, new_filename, new_header, format=False):
        genome_length=0
        with gzip.open(filename, 'rt') as old, gzip.open(new_filename, 'wt') as new:
            for seq_record in SeqIO.parse(old, 'fasta'):
                seq_record.id=new_header
                genome_length += len(seq_record)
                SeqIO.write(seq_record, new, 'fasta')
        if not format:
            os.remove(filename)
        return genome_length


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


def oldies_cleaner(new_genomes, old, oldies_path):
    # When a genome with same accession prefix but different version is downloaded, the old version is deleted
    # from the current database
    old_no_version=list(map(lambda genome: genome.split(sep='.')[0], old))
    genomes_length=pickle.load(open(os.path.join(oldies_path, 'genomes_length.pkl'), 'rb'))
    for genome, genome_no_version in zip(old, old_no_version):
        if genome_no_version in new_genomes:
            os.remove(os.path.join(oldies_path, genome))
            accession = genome[:-7].split(sep='_')[-1]
            genomes_length.pop(accession)
            print(f'Removing {genome}, new version found')
    pickle.dump(genomes_length, open(os.path.join(oldies_path, 'genomes_length.pkl'), 'wb'))
