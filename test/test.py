import time

import test_database, test_fetcher, test_aligner

if __name__=='__main__':

    time0=time.time()
    t=test_fetcher.test_ftp_selector(mode='single', species=['Xylella'])
    print(f'Table processing and ftp selection took {time.time()-time0} seconds')

    time0=time.time()
    oldies=test_fetcher.test_fetcher(t)
    print(f'Genomes download took {time.time()-time0} seconds')

    time0=time.time()
    db=test_database.test_database_builder(oldies)
    print(f'Database formation took {time.time()-time0} seconds')

    time0=time.time()
    idx=test_aligner.test_aligner_indexer(db)
    print(f'Index building took {time.time()-time0} seconds')

    #time0 = time.time()
    #test_aligner.test_aligner_indexer(idx_file=True)
    #print(f'Index building with file took {time.time()-time0} seconds')