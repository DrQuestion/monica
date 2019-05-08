import os
import shutil
import pickle

from .fetcher import GENOMES_PATH


GENOMES = os.path.join(GENOMES_PATH, '*.fna.gz')
DATABASE_PATH = os.path.join(GENOMES_PATH, 'database')
DATABASE_NAME = ['database', '.fna.gz']


def multi_threaded_builder(genomes=None, database_name=DATABASE_NAME, keep_genomes=None, oldies_path=None, n_threads=None):
    if not os.path.exists(DATABASE_PATH):
        os.mkdir(DATABASE_PATH)
    else:
        for database in os.listdir(DATABASE_PATH):
            if database.endswith('.fna.gz'):
                os.remove(database)

    current_genomes_length=dict()

    if keep_genomes:

        # update genomes_length / create ex-novo
        current_genomes_length = pickle.load(open(os.path.join(GENOMES_PATH, 'current_genomes_length.pkl'), 'rb'))
        if 'genomes_length.pkl' in os.listdir(oldies_path):
            genomes_length = pickle.load(open(os.path.join(oldies_path, 'genomes_length.pkl'), 'rb'))
            genomes_length.update(current_genomes_length)
        else:
            genomes_length = current_genomes_length
        pickle.dump(genomes_length, open(os.path.join(oldies_path, 'genomes_length.pkl'), 'wb'))

        # move genomes in oldies_path
        for database in os.listdir(GENOMES_PATH):
            if database.endswith('.fna.gz'):
                shutil.move(os.path.join(GENOMES_PATH, database), oldies_path)
                # would raise errors if same genome downloaded twice, but it should not happen for monica's structure
                # looks for an old genome before downloading it
    else:
        # delete genomes without storing them
        for database in os.listdir(GENOMES_PATH):
            if database.endswith('.fna.gz'):
                os.remove(os.path.join(GENOMES_PATH, database))

    with open(os.path.join(GENOMES_PATH, 'database_created'), 'wb'):
        pass
    return DATABASE_PATH


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
                  .format(genome[0], genome[1][0], (max_chunk_size - size)*16))
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

