from context import fetcher

def test_ftp_selector(mode=None, species=[]):
    t=fetcher.ftp_selector(mode, species)
    return t

def test_fetcher(table=None):
    fetcher.fetcher(table=table)

if __name__=='__main__':
    t=test_ftp_selector(mode='single', species=['Xylella fastidiosa'])
    test_fetcher(t)