from context import aligner

def test_aligner_indexer(idx_file=False):
    idx=aligner.indexer(idx_file=idx_file)
    return idx