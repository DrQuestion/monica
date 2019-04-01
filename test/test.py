import time

import test_database, test_fetcher

if __name__=='__main__':

    time0=time.time()
    t=test_fetcher.test_ftp_selector(mode='all', species=['Xylella fastidiosa'])
    print(f'Table processing and ftp selection took {time.time()-time0} seconds')

    time0=time.time()
    oldies=test_fetcher.test_fetcher(t)
    print(f'Genomes download took {time.time()-time0} seconds')

    time0=time.time()
    test_database.test_database_builder(oldies)
    print(f'Database formation took {time.time()-time0} seconds')

    # time1=time.time()
    # test_database.test_database_builder_alt()
    # print(time.time()-time1)