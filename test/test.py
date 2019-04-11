import os
import time
import shutil

import test_aligner, test_database, test_fetcher

if __name__ == '__main__':

    try:
        os.path.exists('/home/drq/Desktop/temp/query')
    except:
        print('Error')


    query_folder='/home/drq/Desktop/temp/query'

    time0=time.time()
    t=test_fetcher.test_ftp_selector(mode='all', species=['Xylella'])
    print(f'Table processing and ftp selection took {time.time()-time0} seconds')

    time0=time.time()
    oldies=test_fetcher.test_fetcher(t)
    print(f'Genomes download took {time.time()-time0} seconds')

    time0=time.time()
    db=test_database.test_database_builder(oldies)
    print(f'Database formation took {time.time()-time0} seconds')

    # shutil.rmtree(os.path.join(query_folder, 'mapped'))
    # shutil.rmtree(os.path.join(query_folder, 'unmapped'))

    time0=time.time()
    idx=test_aligner.test_aligner_indexer(db, n_threads=7)
    print(f'Index building 3 threads took {time.time()-time0} seconds')

    time0 = time.time()
    print('started alignment')
    alignment = test_aligner.test_aligner_aligner(query_folder, idx)
    print(f'Alignment 3 threads took {time.time()-time0} seconds')

    # time0=time.time()
    # alignment = test_aligner.test_aligner_normalizer(alignment)
    # print(f'Normalization took {time.time()-time0} seconds')
    #
    # time0=time.time()
    # df= test_aligner.test_aligner_alignment_to_data_frame(alignment)
    # print(f'Dataframe conversion took {time.time()-time0} seconds')
