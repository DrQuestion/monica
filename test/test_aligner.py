import pickle
from context import aligner


def test_aligner_indexer(database, idx_file=False, n_threads=None):
    idx = aligner.indexer(database, idx_file=idx_file, n_threads=n_threads)
    return idx


def test_aligner_aligner(query_folder, index, mode='basic'):
    alignment = aligner.aligner(query_folder, index, mode=mode)
    pickle.dump(alignment, open('/home/drq/Desktop/al.pkl', 'wb'))
    return alignment


def test_aligner_normalizer(alignment):
    alignment = aligner.normalizer(alignment)
    pickle.dump(alignment, open('/home/drq/Desktop/norm_al.pkl', 'wb'))
    return alignment


def test_aligner_alignment_to_data_frame(alignment):
    df = aligner.alignment_to_data_frame(alignment)
    df.to_pickle('/home/drq/Desktop/df.pkl')
    return df
