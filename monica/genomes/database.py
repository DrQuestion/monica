import os
import shutil
import pickle

from .fetcher import GENOMES_PATH


GENOMES=os.path.join(GENOMES_PATH, '*.fna.gz')
DATABASE_PATH=os.path.join(GENOMES_PATH, 'database')
DATABASE_NAME='database.fna.gz'


def builder(oldies, database_name=DATABASE_NAME, keep_genomes=None, oldies_path=None):
    database = os.path.join(DATABASE_PATH, database_name)
    if not os.path.exists(DATABASE_PATH):
        os.mkdir(DATABASE_PATH)
        os.system(f'cat {GENOMES} >> {database}')
    else:
        if os.path.isfile(database):
            os.remove(database)
        if oldies:
            for oldie in oldies:
                os.system(f'cat {os.path.join(oldies_path, oldie)} >> {database}')
        if any(fname.endswith('.fna.gz') for fname in os.listdir(GENOMES_PATH)):
            os.system(f'cat {GENOMES} >> {database}')

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
        for file in os.listdir(GENOMES_PATH):
            if file.endswith('.fna.gz'):
                shutil.move(os.path.join(GENOMES_PATH, file), oldies_path)
                # would raise errors if same genome downloaded twice, but it should not happen for monica's structure
                # looks for an old genome before downloading it
    else:
        # delete genomes without storing them
        for file in os.listdir(GENOMES_PATH):
            if file.endswith('.fna.gz'):
                os.remove(os.path.join(GENOMES_PATH, file))

    return database
