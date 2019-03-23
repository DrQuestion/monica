from context import fetcher,tables
import os
import pandas as pd



def test_tables_fetcher():
    try:
        tables.fetcher()
        print('\n----------------------------------\ntest_tables_fetcher worked\n----------------------------------\n')
    except:
        print('test_tables_fetcher broke')


def test_tables_updated():
    try:
        os.chdir(tables.TABLES_PATH)
        ret=tables.updated()
        with open('log', 'r') as log:
            date = log.read()
        os.chdir(tables.CWD)
        print(f'----------------------------------\ntest_tables_updated worked, returned {ret}\nlog content is {date}\n----------------------------------\n')
    except:
        print('test_tables_updated broke')


def test_tables_importer():
    #null case
    try:
        ret=tables.importer()
        print(f'----------------------------------\ntest_tables_importer, which=None, worked, returned {ret}\n----------------------------------\n')
    except:
        print('test_tables_importer, which=None, broke')

    #refseq case
    try:
        ret=tables.importer(which='refseq')
        ret=pd.DataFrame.head(ret)
        print(f'----------------------------------\ntest_tables_importer, which=refseq, worked, returned\n\n{ret}\n----------------------------------\n')
    except:
        print('----------------------------------\ntest_tables_importer, which=refseq, broke\n----------------------------------\n')

    #genbank case
    try:
        ret=tables.importer(which='genbank')
        ret=pd.DataFrame.head(ret)
        print(f'----------------------------------\ntest_tables_importer, which=genbank, worked, returned\n\n{ret}\n----------------------------------\n')
    except:
        print('----------------------------------\ntest_tables_importer, which=genbank, broke\n----------------------------------\n')



if __name__=='__main__':

    test_tables_fetcher()

    test_tables_updated()

    test_tables_importer()
