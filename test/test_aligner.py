import os
from context import aligner
from monica.plots import barplot


def test_aligner_indexer(databases):
    indexes = aligner.indexer(databases=databases)
    return indexes


def test_aligner_multi_threaded_aligner(query_folder, indexes, mode='basic', n_threads=None, focus=None, output=None):
    alignment = aligner.multi_threaded_aligner(query_folder, indexes, mode=mode, n_threads=n_threads,
                                               focus_species=focus, output_folder=output)
    return alignment


def test_aligner_normalizer(alignment, genomes_length=None):
    alignment = aligner.normalizer(alignment, genomes_length=genomes_length)
    return alignment


def test_aligner_alignment_to_data_frame(alignment, output=None):
    df = aligner.alignment_to_data_frame(alignment, output_folder=output)
    return df


if __name__ == '__main__':

    query = '/data/aalbaneseData/monica_data/query_Xylella03/test/test'
    output = '/data/aalbaneseData/monica_data/query_Xylella03/test/output_test'
    if not os.path.exists(output):
        os.mkdir(output)

    # dbs = aligner.DATABASES_PATH
    # idxs = test_aligner_indexer(dbs)
    # indexes_paths = [os.path.join(aligner.INDEXES_PATH, file) for file in
    #                 os.listdir('/data/aalbaneseData/monica_data/query_Xylella03/test/indexes')]
    indexes_paths = [os.path.join(aligner.INDEXES_PATH, file) for file in
                     os.listdir(aligner.INDEXES_PATH) if file.endswith('.mmi')]

    print(indexes_paths)

    alignment = test_aligner_multi_threaded_aligner(query, indexes_paths, mode='query_length', n_threads=6,
                                                    output=output, focus=['Xylella_fastidiosa'])

    print(alignment)

    raw_alignment_df = aligner.alignment_to_data_frame(alignment, output_folder=output,
                                                        filename='raw_monica.dataframe')

    norm_al = test_aligner_normalizer(alignment)

    al_df = test_aligner_alignment_to_data_frame(alignment, output=output)

    barplot.plotter(al_df, raw_alignment_df, output_folder=output, palette='jet', mode='single',
                    show_legend=False, auto_open=False)

    #barplot.plotter(alignment, al_df, output_folder=output, mode='single', palette='jet',
                    #guests='Xanthomonadaceae', show_legend=False)
