import os

from context import database, fetcher


def test_database_multi_threaded_builder(genomes=None, max_chunk_size=None, database_name=None, keep_genomes=None, n_threads=None):
    if database_name:
        db = database.multi_threaded_builder(genomes, max_chunk_size=max_chunk_size, database_name=database_name, keep_genomes=keep_genomes, n_threads=n_threads)
    else:
        db = database.multi_threaded_builder(genomes, max_chunk_size=max_chunk_size, keep_genomes=keep_genomes, n_threads=n_threads)
    return db


def test_database__genomes_slpitter(genomes=None, max_chunk_size=None):
    return database._genomes_splitter(genomes=genomes, max_chunk_size=max_chunk_size)


if __name__ == '__main__':

    MAX = 0.014*(2**30)

    t = fetcher.ftp_selector(mode='single', species=['Xanthomonadaceae'])
    genomes = fetcher.fetcher(table=t, keep_genomes=True)

    for chunk in test_database__genomes_slpitter(genomes=genomes, max_chunk_size=MAX):
        print('A chunk is {}\n'.format(chunk))

    test_database_multi_threaded_builder(genomes=genomes, max_chunk_size=MAX, keep_genomes=True, n_threads=8)
