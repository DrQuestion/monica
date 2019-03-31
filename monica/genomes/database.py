import os
#import gzip
from monica.genomes.fetcher import GENOMES_PATH

DATABASE_PATH=os.path.join(GENOMES_PATH, 'database')
CWD=os.getcwd()
DATABASE_NAME='database.fna.gz'

def builder():
    if not os.path.exists(DATABASE_PATH):
        os.mkdir(DATABASE_PATH)

    os.chdir(GENOMES_PATH)
    database=os.path.join(DATABASE_PATH,DATABASE_NAME)
    os.system(f'cat *.fna.gz >> {database}')
    os.chdir(CWD)

# def builder_alt():
#     os.chdir(GENOMES_PATH)
#     genomes=os.listdir('.')
#     with gzip.open(DATABASE_NAME, 'wt') as database:
#         for genome in genomes:
#             with gzip.open(genome, 'rt') as current_genome:
#                 for line in current_genome:
#                     database.write(line)
#     os.chdir(CWD)

