import os
import datetime as dt
import pandas as pd
import wget
import time

with open(os.path.join(os.path.join(os.path.expanduser('~'), '.monica'), '.root'), 'r') as root:
    MONICA_ROOT = root.readline()

TABLES_PATH = os.path.join(MONICA_ROOT, 'tables')
REFSEQ_SUMMARY_FTP = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt'
GENBANK_SUMMARY_FTP = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt'
LOG_FILE = os.path.join(TABLES_PATH, 'log')
HEADER_LINE = 1
USE_COLS = [0, 5, 6, 7, 8, 9, 19]


def fetcher():
    start = time.time()
    if not os.path.exists(TABLES_PATH):
        os.mkdir(TABLES_PATH)
    if not updated():
        with open(LOG_FILE, 'w') as log:
            log.write(str(dt.date.today()))
        print('Started fetching tables')
        wget.download(REFSEQ_SUMMARY_FTP, out=TABLES_PATH)
        wget.download(GENBANK_SUMMARY_FTP, out=TABLES_PATH)
        print('\nFinished fetching tables')
    print(f'Tables download took {time.time()-start} seconds')


def importer(which=None):
    fetcher()
    if which == 'refseq':
        refseq_table = pd.read_csv(f'{TABLES_PATH}/{REFSEQ_SUMMARY_FTP.split(sep="/")[-1]}', header=HEADER_LINE,
                                   dtype='str', sep='\t', usecols=USE_COLS)
        return refseq_table
    elif which == 'genbank':
        genbank_table = pd.read_csv(f'{TABLES_PATH}/{GENBANK_SUMMARY_FTP.split(sep="/")[-1]}', header=HEADER_LINE,
                                    dtype='str', sep='\t', usecols=USE_COLS)
        return genbank_table
    return 0


def updated():
    if not os.listdir(TABLES_PATH):
        return 0
    with open(LOG_FILE, 'r') as log:
        date = log.read()
    date = dt.datetime.strptime(date, '%Y-%m-%d')
    delta = dt.datetime.now() - date
    if delta.days > 2:
        os.remove(os.path.join(TABLES_PATH, REFSEQ_SUMMARY_FTP.split(sep='/')[-1]))
        os.remove(os.path.join(TABLES_PATH, GENBANK_SUMMARY_FTP.split(sep='/')[-1]))
        return 0
    return 1
