import pickle
from context import aligner


def test_aligner_indexer(database, idx_file=False, n_threads=None):
    idx = aligner.indexer(database, idx_file=idx_file, n_threads=n_threads)
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
