import os
import shutil
from monica.genomes.fetcher import GENOMES_PATH, OLDIES_PATH


GENOMES=os.path.join(GENOMES_PATH, '*.fna.gz')
DATABASE_PATH=os.path.join(GENOMES_PATH, 'database')
DATABASE_NAME='database.fna.gz'


def builder(oldies, database_name=DATABASE_NAME, keep_genomes=True):
    database = os.path.join(DATABASE_PATH, database_name)
    if not os.path.exists(DATABASE_PATH):
        os.mkdir(DATABASE_PATH)
        os.system(f'cat {GENOMES} >> {database}')
    else:
        os.remove(database)
        if oldies:
            for oldie in oldies:
                os.system(f'cat {os.path.join(OLDIES_PATH, oldie)} >> {database}')
        if any(fname.endswith('.fna.gz') for fname in os.listdir(GENOMES_PATH)):
            os.system(f'cat {GENOMES} >> {database}')

    if keep_genomes:
        for file in os.listdir(GENOMES_PATH):
            if file.endswith('.fna.gz'):
                shutil.move(os.path.join(GENOMES_PATH, file), OLDIES_PATH)
                # would raise errors if same genome downloaded twice, but it should not happen for monica's structure
                # looks for an old genome before downloading it
    else:
        for file in os.listdir(GENOMES_PATH):
            if file.endswith('.fna.gz'):
                os.remove(file)

    return database
