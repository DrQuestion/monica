from context import fetcher

def test_ftp_selector(mode=None, species=[]):
    t=fetcher.ftp_selector(mode, species)
    return t

def test_fetcher(table=None):
    oldies=fetcher.fetcher(table=table)
    return oldies

def test_oldies_cleaner():
    fetcher.oldies_cleaner()

if __name__=='__main__':
    test_oldies_cleaner()