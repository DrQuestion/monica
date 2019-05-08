from context import fetcher


def test_ftp_selector(mode=None, species=[]):
    t = fetcher.ftp_selector(mode, species)
    return t


def test_fetcher(table=None, oldies_path=fetcher.OLDIES_PATH, keep_genomes=None, format_genomes=None):
    genomes = fetcher.fetcher(table=table, oldies_path=oldies_path,
                              keep_genomes=keep_genomes, format_genomes=format_genomes)
    return genomes


def main_test ():
    t = test_ftp_selector(mode='single', species=['Xanthomonadaceae'])
    genomes = test_fetcher(table=t, keep_genomes=True)


if __name__=='__main__':
    main_test()
