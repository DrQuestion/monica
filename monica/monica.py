#! /home/drq/.venvs/monica/bin/python

import os
import argparse
import pickle

import psutil

import genomes.fetcher as gfetcher
import genomes.database as gdatabase
import genomes.aligner as galigner
import plots.barplot as barplot
import helpers.helpers as helpers


def main():
    parser = argparse.ArgumentParser(prog='monica')

    input_folder = parser.add_mutually_exclusive_group(required=True)
    input_folder.add_argument('-q', '--query_folder',
                              help='The folder where input fastq are stored')
    input_folder.add_argument('-f5', '--fast5_folder',
                              help='The folder where input fast5 are stored (real-time use only)')
    parser.add_argument('-o', '--output_folder',
                        help='The folder where output files will be stored')
    parser.add_argument('-t', '--threads', type=int, default=3,
                        help='The number of threads to be used')

    subparsers = parser.add_subparsers()

    from_scratch = subparsers.add_parser('from_scratch', help='To run monica from the beginning')
    from_scratch.add_argument('-H', '--host_specie',
                        help='Host species where samples come from, underscore separated ("Genus_specie")')
    from_scratch.add_argument('-G', '--guest_species', nargs='*',
                        help='Guest species hypothesized to be present in the sample, '
                             'underscore separated ("Genus_specie"), but also higher levels of the tree of life, '
                             'like order or genus only')
    from_scratch.add_argument('-k', '--keep_genomes', choices=['yes', 'no'], required=True,
                        help='Choose if genomes are to be kept in -g')
    from_scratch.add_argument('-g', '--genomes_folder', default=gfetcher.OLDIES_PATH,
                        help='The folder where already stored genomes will be found '
                             'and new genomes will be stored if --keep_genomes option is set to "yes"')
    from_scratch.add_argument('--format_genomes',
                        help='Folder, call it together with -m and -H/-G, when assembly-formatted genomes '
                             'are present in it without being formatted for monica')
    from_scratch.add_argument('-im', '--indexing_memory', type=helpers.human_readable,
                        help='The maximum memory the user wishes to be used during indexing. '
                             'Must end with a unit of measure among B|K|M|G|T.'
                             'Default is half of total memory')
    from_scratch.add_argument('-m', '--mode', choices=['single', 'all', 'overnight'],
                              help="The mode of monica's execution: single, all or overnight")
    from_scratch.set_defaults(func=main_from_scratch)

    from_alignment = subparsers.add_parser('from_alignment', help='To run monica given an already built index/database')
    from_alignment.add_argument('-uoi', '--use_old_index', action='store_true',
                        help='Use the index of the previous analysis, useful for frequent use of the same genomes')
    # momentary option, active since -uoi will be stable:
    from_alignment.add_argument('-uod', '--use_old_database', action='store_true')
    from_alignment.set_defaults(func=main_from_alignment)

    args = parser.parse_args()
    args.func(args)
    return args


def main_from_scratch(args):

    # Input handling
    if args.query_folder:
        input_folder = args.query_folder
    else:
        input_folder = args.fast5_folder
    if args.output_folder:
        output_folder = args.output_folder
    else:
        output_folder = os.path.join(input_folder, 'monica_output')
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    if args.indexing_memory:
        indexing_memory = args.indexing_memory
    else:
        indexing_memory = psutil.virtual_memory().total/2
    n_threads = args.threads

    mode = args.mode
    host = args.host_specie
    guests = args.guest_species
    if guests:
        guests = map(lambda guest: ' '.join(guest.split(sep='_')), guests)
    if args.keep_genomes == 'yes':
        keep_genomes = True
    else:
        keep_genomes = False
    oldies_path = args.genomes_folder
    format_genomes = args.format_genomes

    # Genomes ftp selection and download
    ftp_table = gfetcher.ftp_selector(mode=mode, species=guests)

    if host:
        host = ' '.join(host.split(sep='_'))
        host_ftp_table = gfetcher.ftp_selector(mode='single', species=[host])
        ftp_table = ftp_table.append(host_ftp_table, ignore_index=True)

    oldies = gfetcher.fetcher(ftp_table, oldies_path=oldies_path, format_genomes=format_genomes)

    # Database building and indexing
    database = gdatabase.builder(oldies, keep_genomes=keep_genomes, oldies_path=oldies_path)

    with open(os.path.join(gfetcher.GENOMES_PATH, 'going_to_enter_indexing'), 'wb'):
        pass

    idx = galigner.indexer(database, n_threads=n_threads)

    with open(os.path.join(gfetcher.GENOMES_PATH, 'going_to_enter_alignment'), 'wb'):
        pass

    # Alignment and normalization
    alignment = galigner.multi_threaded_aligner(input_folder, idx, mode='basic', n_threads=n_threads)

    norm_alignment = galigner.normalizer(alignment)

    alignment_df = galigner.alignment_to_data_frame(norm_alignment, output_folder=output_folder)

    # Plotting
    barplot.plotter(norm_alignment, alignment_df, host=host, output_folder=output_folder)


def main_from_alignment(args):

    # Input handling
    if args.query_folder:
        input_folder = args.query_folder
    else:
        input_folder = args.fast5_folder
    if args.output_folder:
        output_folder = args.output_folder
    else:
        output_folder = os.path.join(input_folder, 'monica_output')
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    n_threads = args.threads

    use_old_index = args.use_old_index
    use_old_database = args.use_old_database

    if use_old_database:
        database = os.path.join(gdatabase.DATABASE_PATH, gdatabase.DATABASE_NAME)

        with open(os.path.join(gfetcher.GENOMES_PATH, 'going_to_enter_indexing'), 'wb'):
            pass

        idx = galigner.indexer(database, n_threads=n_threads)

    elif use_old_index:
        idx = pickle.load(open(galigner.INDEX_PICKLE, 'r'))

    with open(os.path.join(gfetcher.GENOMES_PATH, 'going_to_enter_alignment'), 'wb'):
        pass

    # Alignment and normalization
    alignment = galigner.multi_threaded_aligner(input_folder, idx, mode='basic', n_threads=n_threads)

    norm_alignment = galigner.normalizer(alignment)

    alignment_df = galigner.alignment_to_data_frame(norm_alignment, output_folder=output_folder)

    # Plotting
    barplot.plotter(norm_alignment, alignment_df, output_folder=output_folder)


if __name__ == '__main__':

    main()


