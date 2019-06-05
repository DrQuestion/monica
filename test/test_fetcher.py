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


if __name__ == '__main__':

    t = test_ftp_selector(mode='focus', species=['Xylella fastidiosa'])
    fetcher.focus_fetcher(t)

# ['2371', '1947786', '1444770', '405952', '405953', '405954', '570946', '542288', '483031', '380250']
# ['2371', '1947786', '1444770', '405952', '405953', '405954', '570946', '542288', '483031', '380250']
# ['155919', '160492', '183190', '329530', '405440', '405441', '788929', '1123506', '1395102', '1267006', '1343737',
# '155920', '1445511', '945689', '1211847', '1214121', '1401256', '1403344', '1453499', '380250', '405952', '405953',
# '405954', '483031', '542288', '570946', '1444770', '1947786']
