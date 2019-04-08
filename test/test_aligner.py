from context import aligner

def test_aligner_indexer(database, idx_file=False):
    idx=aligner.indexer(database, idx_file=idx_file)
    return idx