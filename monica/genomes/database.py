import os
import shutil
from monica.genomes.fetcher import GENOMES_PATH, OLDIES_PATH


GENOMES=os.path.join(GENOMES_PATH, '*.fna.gz')
DATABASE_PATH=os.path.join(GENOMES_PATH, 'database')
DATABASE_NAME='database.fna.gz'
DATABASE=os.path.join(DATABASE_PATH,DATABASE_NAME)


def builder(oldies, keep_genomes=True):
    if not os.path.exists(DATABASE_PATH):
        os.mkdir(DATABASE_PATH)
        os.system(f'cat {GENOMES} >> {DATABASE}')
    else:
        os.remove(DATABASE)
        if oldies:
            for oldie in oldies:
                os.system(f'cat {os.path.join(OLDIES_PATH, oldie)} >> {DATABASE}')
        if any(fname.endswith('.fna.gz') for fname in os.listdir(GENOMES_PATH)):
            os.system(f'cat {GENOMES} >> {DATABASE}')

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
