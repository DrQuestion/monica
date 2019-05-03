import os
from context import aligner, database


def test_aligner_indexer(database, n_threads=None):
    idx = aligner.indexer(database, n_threads=n_threads)
    return idx


def test_aligner_multi_threaded_aligner(query_folder, index, mode='basic', n_threads=None):
    alignment = aligner.multi_threaded_aligner(query_folder, index, mode=mode, n_threads=n_threads)
    return alignment


def test_aligner_normalizer(alignment):
    alignment = aligner.normalizer(alignment)
    return alignment


def test_aligner_alignment_to_data_frame(alignment):
    df = aligner.alignment_to_data_frame(alignment)
    return df


if __name__ == '__main__':

    db = os.path.join(aligner.DATABASE_PATH, database.DATABASE_NAME)
    idx = test_aligner_indexer(db, n_threads=4)
