import time

import test_aligner, test_database, test_fetcher
from monica.plots import barplot

if __name__ == '__main__':

    max_chunk_size = 0.015 * (2 ** 30)
    query_folder = '/home/drq/temp_query'
    n_threads = 8
    output = '/home/drq/temp_output'
    guests = ['Xanthomonadaceae', 'Aquimonas', 'Dokdonella']
    mode = 'single'
    show_legend = False
    auto_open = True

    time0 = time.time()
    t = test_fetcher.test_ftp_selector(mode=mode, species=guests)
    print(f'Table processing and ftp selection took {time.time()-time0} seconds')

    time1 = time.time()
    genomes = test_fetcher.test_fetcher(table=t, keep_genomes=True)
    print(f'Genomes download took {time.time()-time1} seconds')

    time2 = time.time()
    dbs, lengths = test_database.test_database_multi_threaded_builder(genomes=genomes, keep_genomes=True,
                                                                      max_chunk_size=max_chunk_size, n_threads=n_threads)
    print(f'Database formation took {time.time()-time2} seconds')
    print(lengths)

    time3 = time.time()
    idxs = test_aligner.test_aligner_indexer(dbs)
    print(f'Index building took {time.time()-time3} seconds')

    time4 = time.time()
    alignment = test_aligner.test_aligner_multi_threaded_aligner(query_folder, idxs, mode='basic', n_threads=n_threads)
    print(f'Alignment took {time.time()-time4} seconds')

    time5 = time.time()
    # alignment=pickle.load(open('/home/drq/PycharmProjects/tesi/monica/monica/genomes/alignment.pkl', 'rb'))
    alignment = test_aligner.test_aligner_normalizer(alignment, lengths)
    print(f'Normalization took {time.time()-time5} seconds')

    time6 = time.time()
    df= test_aligner.test_aligner_alignment_to_data_frame(alignment, output=output)
    print(f'Dataframe conversion took {time.time()-time6} seconds')

    barplot.plotter(norm_alignment=alignment, norm_alignment_df=df, output_folder=output, guests=guests, mode=mode,
                    auto_open=auto_open, show_legend=show_legend)

