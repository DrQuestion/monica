import os
import pickle
import gzip
from itertools import count, repeat
from multiprocessing.dummy import Pool as ThreadPool

from Bio import SeqIO

from .fetcher import GENOMES_PATH


DATABASES_PATH = os.path.join(GENOMES_PATH, 'databases')
DATABASE_NAME = ['database', '.fna.gz']


def multi_threaded_builder(genomes=None, max_chunk_size=None, databases_path=DATABASES_PATH,
                           database_name=DATABASE_NAME, keep_genomes=None, n_threads=None):
    if not os.path.exists(databases_path):
        os.mkdir(databases_path)
    else:
        for database in os.listdir(databases_path):
            if database.endswith('.fna.gz'):
                os.remove(os.path.join(databases_path, database))

    current_genomes_length = dict()

    pool = ThreadPool(n_threads)

    lengths = pool.starmap(builder, zip(_genomes_splitter(genomes, max_chunk_size=max_chunk_size),
                                        repeat(databases_path), repeat(database_name), count()))

    for length in lengths:
        current_genomes_length.update(length)

    if not keep_genomes:
        # delete genomes without storing them
        # TODO: if focus and a genome belongs to focus, don't delete it
        for genome in os.listdir(GENOMES_PATH):
            if genome.endswith('.fna.gz'):
                os.remove(os.path.join(GENOMES_PATH, genome))

    pickle.dump(current_genomes_length, open(os.path.join(databases_path, 'current_genomes_length.pkl'), 'wb'))

    with open(os.path.join(GENOMES_PATH, 'database_created'), 'wb'):
        pass
    return databases_path, current_genomes_length


def builder(genomes_chunk, databases_path, database_name, database_number):
    database_file = os.path.join(databases_path, str(database_number).join(database_name))
    print('Working on {}'.format(str(database_number).join(database_name)))
    this_database_genomes_length = dict()
    with gzip.open(database_file, 'wt') as database:
        for genome in genomes_chunk:
            genome_length = 0
            new_header = ':'.join(genome[1])
            with gzip.open(genome[0], 'rt') as g:
                for seq_record in SeqIO.parse(g, 'fasta'):
                    seq_record.id = new_header
                    genome_length += len(seq_record)
                    SeqIO.write(seq_record, database, 'fasta')
            this_database_genomes_length[genome[1][1]] = genome_length
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
