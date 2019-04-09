import sys
import getopt

import monica.genomes.fetcher as gfetcher
import monica.genomes.database as gdatabase
import monica.genomes.aligner as galigner


def main(argv):
    opts, args = getopt.getopt(argv[1:], '')
    pass


if __name__=='__main__':
    main(sys.argv)
