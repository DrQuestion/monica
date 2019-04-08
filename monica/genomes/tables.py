import os
import datetime as dt
import pandas as pd
import wget
import time


TABLES_PATH=os.path.join(os.path.dirname(__file__), 'tables')
REFSEQ_SUMMARY_FTP= 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt'
GENBANK_SUMMARY_FTP= 'ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt'
HEADERLINE=1
USECOLS=[0,5,7,8,9,19]
CWD=os.getcwd()


def fetcher():
    start=time.time()
    if not os.path.exists(TABLES_PATH):
        os.mkdir(TABLES_PATH)
    os.chdir(TABLES_PATH)
    if not updated():
        with open('log', 'w+') as log:
            log.write(str(dt.date.today()))
        print('started_fetching_tables')
        wget.download(REFSEQ_SUMMARY_FTP)
        wget.download(GENBANK_SUMMARY_FTP)
        print('finished_fetching_tables')
    os.chdir(CWD)
    print(f'Tables download took {time.time()-start} seconds')

def importer(which=None):
    fetcher()
    if which=='refseq':
        refseq_table = pd.read_csv(f'{TABLES_PATH}/{REFSEQ_SUMMARY_FTP.split(sep="/")[-1]}', header=HEADERLINE, dtype='str', sep='\t', usecols=USECOLS)
        return refseq_table
    elif which=='genbank':
        genbank_table = pd.read_csv(f'{TABLES_PATH}/{GENBANK_SUMMARY_FTP.split(sep="/")[-1]}', header=HEADERLINE, dtype='str', sep='\t', usecols=USECOLS)
        return genbank_table
    return 0


def updated():
    if not os.listdir(TABLES_PATH):
        return 0
    with open('log', 'r') as log:
        date = log.read()
    date = dt.datetime.strptime(date, '%Y-%m-%d')
    delta = dt.datetime.now() - date
    if delta.days > 2:
        os.remove(REFSEQ_SUMMARY_FTP.split(sep='/')[-1])
        os.remove(GENBANK_SUMMARY_FTP.split(sep='/')[-1])
        return 0
    return 1
