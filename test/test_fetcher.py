from context import fetcher
def test_ftp_selector(mode=None, species=[]):
    fetcher.ftp_selector(mode, species)

if __name__=='__main__':
    test_ftp_selector(mode='all', species=['Xylella'])