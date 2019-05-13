import os
from context import aligner
from monica.plots import barplot


def test_aligner_indexer(databases):
    indexes = aligner.indexer(databases=databases)
    return indexes


def test_aligner_multi_threaded_aligner(query_folder, indexes, mode='basic', n_threads=None):
    alignment = aligner.multi_threaded_aligner(query_folder, indexes, mode=mode, n_threads=n_threads)
    return alignment


def test_aligner_normalizer(alignment, genomes_length=None):
    alignment = aligner.normalizer(alignment, genomes_length=genomes_length)
    return alignment


def test_aligner_alignment_to_data_frame(alignment, output=None):
    df = aligner.alignment_to_data_frame(alignment, output_folder=output)
    return df


if __name__ == '__main__':

    query = '/home/drq/temp_query'
    output = '/home/drq/temp_output'
    if not os.path.exists(output):
        os.mkdir(output)

    dbs = aligner.DATABASES_PATH
    # idxs = test_aligner_indexer(dbs)
    idxs = aligner.indexes_opener(aligner.INDEXES_PATH)

    for idx in idxs:
        print(idx)

    alignment = test_aligner_multi_threaded_aligner(query, idxs, mode='basic', n_threads=6)

    print(alignment)

    norm_al = test_aligner_normalizer(alignment)

    al_df = test_aligner_alignment_to_data_frame(alignment, output=output)

    barplot.plotter(alignment, al_df, output_folder=output)
