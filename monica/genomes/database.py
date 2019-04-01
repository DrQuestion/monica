import os
import shutil
#import gzip
from monica.genomes.fetcher import GENOMES_PATH

GENOMES=os.path.join(GENOMES_PATH, '*.fna.gz')
DATABASE_PATH=os.path.join(GENOMES_PATH, 'database')
DATABASE_NAME='database.fna.gz'
DATABASE=os.path.join(DATABASE_PATH,DATABASE_NAME)

def builder(oldies):
    if not os.path.exists(DATABASE_PATH):
        os.mkdir(DATABASE_PATH)
        os.system(f'cat {GENOMES} >> {DATABASE}')

    else:
        os.remove(DATABASE)
        if oldies:
            oldies_path=os.path.join(GENOMES_PATH,'oldies')
            for oldy in oldies:
                       os.system(f'cat {os.path.join(oldies_path, oldy)} >> {DATABASE}')
        if any(fname.endswith('.fna.gz') for fname in os.listdir(GENOMES_PATH)):
            os.system(f'cat {GENOMES} >> {DATABASE}')

    for file in os.listdir(GENOMES_PATH):
        if file.endswith('.fna.gz'):
            shutil.move(os.path.join(GENOMES_PATH, file), os.path.join(GENOMES_PATH, 'oldies'))