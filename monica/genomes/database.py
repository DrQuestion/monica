import os
import pickle
import gzip
import shutil
from itertools import count, repeat
from multiprocessing.dummy import Pool as ThreadPool

from Bio import SeqIO

from .fetcher import GENOMES_PATH


GENOMES = os.path.join(GENOMES_PATH, '*.fna.gz')
DATABASES_PATH = os.path.join(GENOMES_PATH, 'databases')
DATABASE_NAME = ['database', '.fna.gz']
EXCEPTIONS_PATH = os.path.join(GENOMES_PATH, 'exceptions')


def multi_threaded_builder(genomes=None, max_chunk_size=None, database_name=DATABASE_NAME, keep_genomes=None, n_threads=None):
    if not os.path.exists(DATABASES_PATH):
        os.mkdir(DATABASES_PATH)
    else:
        for database in os.listdir(DATABASES_PATH):
            if database.endswith('.fna.gz'):
                os.remove(os.path.join(DATABASES_PATH, database))
    if not os.path.exists(EXCEPTIONS_PATH):
        os.mkdir(EXCEPTIONS_PATH)

    current_genomes_length = dict()

    pool = ThreadPool(n_threads)

    lengths = pool.starmap(builder, zip(_genomes_splitter(genomes, max_chunk_size=max_chunk_size), repeat(database_name), count()))

    for length in lengths:
        current_genomes_length.update(length)

    if not keep_genomes:
        # delete genomes without storing them
        for genome in os.listdir(GENOMES_PATH):
            if genome.endswith('.fna.gz'):
                os.remove(os.path.join(GENOMES_PATH, genome))

    pickle.dump(current_genomes_length, open(os.path.join(GENOMES_PATH, 'current_genomes_length.pkl'), 'wb'))

    with open(os.path.join(GENOMES_PATH, 'database_created'), 'wb'):
        pass
    return DATABASES_PATH, current_genomes_length


def builder(genomes_chunk, database_name, database_number):
    database_file = os.path.join(DATABASES_PATH, str(database_number).join(database_name))
    print('Working on {}'.format(str(database_number).join(database_name)))
    this_database_genomes_length = dict()
    with gzip.open(database_file, 'wt') as database:
        for genome in genomes_chunk:
            genome_length = 0
            new_header = ':'.join(genome[1])
            try:
                with gzip.open(genome[0], 'rt') as g:
                    for seq_record in SeqIO.parse(g, 'fasta'):
                        seq_record.id = new_header
                        genome_length += len(seq_record)
                        SeqIO.write(seq_record, database, 'fasta')
                this_database_genomes_length[genome[1][1]] = genome_length
            except:
                print('{} failed database insertion'.format(genome))
                shutil.move(genome, EXCEPTIONS_PATH)
    return this_database_genomes_length


def _genomes_splitter(genomes, max_chunk_size=None):
    chunk = []
    exceeding_chunk = []
    chunk_size = 0
    for genome in genomes:
        size = os.path.getsize(genome[0])
        if size > max_chunk_size:
            exceeding_chunk.append(genome)
            print('Genome {}, ({}) alone expected to generate an index '
                  'exceeding the maximum memory deriving from settings of {} bytes'
                  .format(genome[0], genome[1][0], (size - max_chunk_size)*16))
            yield exceeding_chunk
            exceeding_chunk = []
        else:
            if chunk_size + size <= max_chunk_size:
                chunk.append(genome)
                chunk_size += size
            else:
                yield chunk
                chunk = []
                chunk_size = 0
    if chunk:
        yield chunk
