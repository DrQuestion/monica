import pickle
from context import aligner


def test_aligner_indexer(database, idx_file=False):
    idx = aligner.indexer(database, idx_file=idx_file)
    return idx


def test_aligner_aligner(query_folder, index, mode='basic'):
    alignment, genomes_length = aligner.aligner(query_folder, index, mode=mode)
    pickle.dump(alignment, open('/home/drq/Desktop/al.pkl', 'wb'))
    pickle.dump(genomes_length, open('/home/drq/Desktop/gl.pkl', 'wb'))
    return alignment, genomes_length


def test_aligner_normalize(alignment, genomes_length):
    alignment = aligner.normalize(alignment, genomes_length)
    pickle.dump(alignment, open('/home/drq/Desktop/norm_al.pkl', 'wb'))
    return alignment


def test_aligner_alignment_to_data_frame(alignment):
    df = aligner.alignment_to_data_frame(alignment)
    df.to_pickle('/home/drq/Desktop/df.pkl')
    return df
