#! /usr/bin/env python3.6

# TODO log system when launching monica from scratch to save all parameters, to keep host, guests and mode
# TODO argcomplete

import os

if os.path.exists(os.path.join(os.path.expanduser('~'), '.monica')):
    with open(os.path.join(os.path.join(os.path.expanduser('~'), '.monica'), '.root'), 'r') as root:
        MONICA_ROOT = root.readline()
        print(f'root is {MONICA_ROOT}')
else:
    os.mkdir(os.path.join(os.path.expanduser('~'), '.monica'))
    with open(os.path.join(os.path.join(os.path.expanduser('~'), '.monica'), '.root'), 'w') as root:
        root.write(os.path.join(os.path.expanduser('~'), '.monica'))
    MONICA_ROOT = os.path.join(os.path.expanduser('~'), '.monica')
os.chdir(MONICA_ROOT)

# import warnings
# warnings.filterwarnings("ignore", category=UserWarning)

import argparse
import psutil
import shutil
import pandas as pd

import genomes.fetcher as gfetcher
import genomes.database as gdatabase
import genomes.aligner as galigner
import plots.barplot as barplot
import helpers.helpers as helpers


def main():
    parser = argparse.ArgumentParser(prog='monica')

    input_folder = parser.add_mutually_exclusive_group()
    input_folder.add_argument('-q', '--query_folder',
                              help='The folder where input fastq are stored')
    input_folder.add_argument('-f5', '--fast5_folder',
                              help='The folder where input fast5 are stored (real-time use only)')
    parser.add_argument('-o', '--output_folder',
                        help='The folder where output files will be stored')
    parser.add_argument('-G', '--guest_species', nargs='*',
                        help='Guest species hypothesized to be present in the sample, '
                             'underscore separated ("Genus_specie"), but also higher levels of the tree of life, '
                             'like order or genus only')
    parser.add_argument('-H', '--host_species', nargs='*',
                        help='Host species where samples come from, underscore separated ("Genus_specie")')
    parser.add_argument('-F', '--focus_species', nargs='*',
                        help='Species of interest among the whole guests selected, to focus analysis on their subspecies')
    parser.add_argument('-m', '--mode', choices=['single', 'all', 'overnight'],
                        help="The mode of monica's execution: single, all or overnight")
    parser.add_argument('-a', '--alignment_mode', default='query_length',
                        help='The alignment mode to be performed.'
                             'Default is %(default)s')
    parser.add_argument('-g', '--genomes_folder', default=gfetcher.OLDIES_PATH,
                        help='The folder where already stored genomes will be found '
                             'and new genomes will be stored if --keep_genomes option is set to "yes". '
                             'Default is %(default)s')
    parser.add_argument('-k', '--keep_genomes', choices=['yes', 'no'], default='yes',
                        help='Choose if genomes are to be kept in -g. Default is "%(default)s"')
    parser.add_argument('-i', '--indexes', nargs='*', help='folder where indexes must be written/read')
    parser.add_argument('--not_auto_open_plot', action='store_true',
                        help='If wishing NOT to plot the results at the and of the analysis')
    parser.add_argument('--not_show_legend', action='store_true',
                        help='If wishing NOT to show legend')
    parser.add_argument('-R', '--reads_threshold', default=0,
                        help='Threshold of reads for a specie in each sample that will be displayed in the plot. '
                             'If lower in all samples species not displayed. '
                             'Default is %(default)d')
    parser.add_argument('--format_genomes',
                        help='Folder, call it together with -m and -H/-G, when assembly-formatted genomes '
                             'are present in it without being formatted for monica')
    parser.add_argument('-t', '--threads', type=int, default=3,
                        help='The number of threads to be used.'
                             'Default is %(default)d')
    parser.add_argument('-im', '--indexing_memory', type=helpers.human_readable,
                        help='The maximum memory the user wishes to be used during indexing. '
                             'Must end with a unit of measure among B|K|M|G|T.'
                             'Default is fourth of total memory')

    subparsers = parser.add_subparsers()

    # from_scratch = subparsers.add_parser('from_scratch', help='To run monica from the beginning')
    # from_scratch.add_argument('-H', '--host_specie',
    #                           help='Host species where samples come from, underscore separated ("Genus_specie")')
    # from_scratch.add_argument('-G', '--guest_species', nargs='*',
    #                           help='Guest species hypothesized to be present in the sample, '
    #                                'underscore separated ("Genus_specie"), but also higher levels of the tree of life, '
    #                                'like order or genus only')
    # from_scratch.add_argument('-F', '--focus_species', nargs='*',
    #                           help='Species of interest among the whole guests selected, focus analysis on their subspecies')
    # from_scratch.add_argument('--format_genomes',
    #                           help='Folder, call it together with -m and -H/-G, when assembly-formatted genomes '
    #                                'are present in it without being formatted for monica')
    # from_scratch.add_argument('-m', '--mode', choices=['single', 'all', 'overnight'],
    #                           help="The mode of monica's execution: single, all or overnight")
    # from_scratch.set_defaults(func=main_from_scratch)

    build_index = subparsers.add_parser('build_index', help='To build an index a priori')
    build_index.set_defaults(func=main_build_index)

    list_indexes = subparsers.add_parser('list_indexes', aliases=['ls', 'li'],
                                         help='To print on the screen all the available built indexes')
    list_indexes.set_defaults(func=main_list_indexes)

    plot_only = subparsers.add_parser('plot_only', help='To generate the plot only, given a table of monica results')
    plot_only.add_argument('-d', '--data_frame', help='The monica result data frame to generate the plot from')
    plot_only.set_defaults(func=main_from_plotting)

    initialize = subparsers.add_parser('initialize', help='To control monica root folder where heavy files will be '
                                                          'downloaded')
    initialize.add_argument('-r', '--monica_root_folder', default=MONICA_ROOT,
                            help='Path where ".monica" root folder will be created')
    initialize.set_defaults(func=main_initialize)

    parser.set_defaults(func=main_after_seq)
    args = parser.parse_args()
    args.func(args)

    return args


# def main_from_scratch(args):
#     # Input handling:
#     if args.query_folder:
#         input_folder = args.query_folder
#     else:
#         input_folder = args.fast5_folder
#     if args.output_folder:
#         output_folder = args.output_folder
#     else:
#         output_folder = os.path.join(input_folder, 'monica_output')
#     if not os.path.exists(output_folder):
#         os.mkdir(output_folder)
#
#     with open(os.path.join(MONICA_ROOT, 'monica.params'), 'w') as params:
#         params.write(str(args))
#
#     n_threads = args.threads
#     alignment_mode = args.alignment_mode
#
#     if not args.focus_species:
#         focus_species = []
#     else:
#         focus_species = args.focus_species
#
#     auto_open_plot = not args.not_auto_open_plot
#     show_legend = not args.not_show_legend
#     reads_threshold = args.reads_threshold
#
#     if args.indexing_memory:
#         indexing_memory = args.indexing_memory
#     else:
#         indexing_memory = psutil.virtual_memory().total / 4
#
#     max_chunk_size = indexing_memory / 16
#
#     mode = args.mode
#     host = args.host_specie
#     guests = args.guest_species
#
#     helpers.initializer(output_folder)
#
#     if guests:
#         guests = map(lambda guest: ' '.join(guest.split(sep='_')), guests)
#     if args.keep_genomes == 'yes':
#         keep_genomes = True
#     else:
#         keep_genomes = False
#     oldies_path = args.genomes_folder
#     format_genomes = args.format_genomes
#
#     # Genomes ftp selection and download:
#     ftp_table = gfetcher.ftp_selector(mode=mode, species=guests)
#
#     if host:
#         host = ' '.join(host.split(sep='_'))
#         host_ftp_table = gfetcher.ftp_selector(mode='single', species=[host])
#         ftp_table = ftp_table.append(host_ftp_table, ignore_index=True)
#
#     genomes = gfetcher.fetcher(ftp_table, oldies_path=oldies_path, keep_genomes=keep_genomes,
#                                format_genomes=format_genomes)
#
#     # Database building and indexing:
#     databases, genomes_length = gdatabase.multi_threaded_builder(genomes=genomes, max_chunk_size=max_chunk_size,
#                                                                  databases_path=gdatabase.DATABASES_PATH,
#                                                                  keep_genomes=keep_genomes, n_threads=n_threads)
#
#     indexes_paths = galigner.indexer(databases)
#
#     # Focus mode indexing building block:
#     if focus_species:
#         focus_input_folder = os.path.join(input_folder, 'focus')
#         focus_output_folder = os.path.join(output_folder, 'focus')
#         if not os.path.exists(focus_output_folder):
#             os.mkdir(focus_output_folder)
#
#         helpers.initializer(focus_output_folder)
#
#         focus_databases_path = os.path.join(gdatabase.DATABASES_PATH, 'focus')
#         focus_species = map(lambda specie: ' '.join(specie.split(sep='_')), focus_species)
#         focus_ftp_table = gfetcher.ftp_selector(mode='focus', species=focus_species)
#         focus_genomes = gfetcher.focus_fetcher(focus_ftp_table, oldies_path=oldies_path, keep_genomes=keep_genomes)
#         focus_databases, focus_genomes_length = gdatabase.multi_threaded_builder(genomes=focus_genomes,
#                                                                                  max_chunk_size=max_chunk_size,
#                                                                                  keep_genomes=keep_genomes,
#                                                                                  n_threads=n_threads,
#                                                                                  databases_path=focus_databases_path)
#         focus_indexes_paths = galigner.indexer(focus_databases,
#                                                indexes_path=os.path.join(galigner.INDEXES_PATH, 'focus'))
#
#     with open(os.path.join(gfetcher.GENOMES_PATH, 'going_to_enter_alignment'), 'wb'):
#         pass
#
#     # Alignment and normalization:
#     alignment = galigner.multi_threaded_aligner(input_folder, indexes_paths, mode=alignment_mode, n_threads=n_threads,
#                                                 focus_species=focus_species, output_folder=output_folder)
#
#     raw_alignment_df = galigner.alignment_to_data_frame(alignment, output_folder=output_folder,
#                                                         filename='raw_monica.dataframe')
#
#     norm_alignment = galigner.normalizer(alignment, genomes_length)
#
#     norm_alignment_df = galigner.alignment_to_data_frame(norm_alignment, output_folder=output_folder)
#
#     # Plotting:
#     barplot.plotter(norm_alignment_df, raw_alignment_df, output_folder=output_folder, palette='jet',
#                     reads_threshold=reads_threshold, hosts=host, guests=args.guest_species, mode=mode,
#                     show_legend=show_legend, auto_open=auto_open_plot)
#
#     if focus_species:
#         focus_alignment = galigner.multi_threaded_aligner(focus_input_folder, focus_indexes_paths, mode=alignment_mode,
#                                                           n_threads=n_threads, output_folder=focus_output_folder)
#         focus_raw_alignment_df = galigner.alignment_to_data_frame(focus_alignment, output_folder=focus_output_folder,
#                                                                   filename='raw_monica.dataframe')
#         focus_norm_alignment = galigner.normalizer(focus_alignment, focus_genomes_length)
#         focus_norm_alignment_df = galigner.alignment_to_data_frame(focus_norm_alignment,
#                                                                    output_folder=focus_output_folder)
#         barplot.plotter(focus_norm_alignment_df, focus_raw_alignment_df, output_folder=focus_output_folder,
#                         palette='jet', reads_threshold=0, guests=focus_species, mode='focus',
#                         show_legend=show_legend, auto_open=auto_open_plot)
#

def main_after_seq(args):
    # Input handling:
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

    with open(os.path.join(MONICA_ROOT, 'monica.params'), 'w') as params:
        params.write(str(args))

    if not args.focus_species:
        focus_species = []
    else:
        focus_species = args.focus_species

    n_threads = args.threads
    alignment_mode = args.alignment_mode

    if args.indexing_memory:
        indexing_memory = args.indexing_memory
    else:
        indexing_memory = psutil.virtual_memory().total / 4

    max_chunk_size = indexing_memory / 16

    if args.keep_genomes == 'yes':
        keep_genomes = True
    else:
        keep_genomes = False
    oldies_path = args.genomes_folder
    format_genomes = args.format_genomes

    auto_open_plot = not args.not_auto_open_plot
    show_legend = not args.not_show_legend
    reads_threshold = args.reads_threshold

    mode = args.mode
    hosts = args.host_species
    guests = args.guest_species

    helpers.initializer(output_folder)

    if guests:
        guest_species = guests
    else:
        guest_species = []
    if hosts:
        host_species = hosts
    else:
        host_species = []

    indexes = args.indexes

    indexes_paths = []
    focus_indexes_paths = []
    inside_indexes = []
    outside_indexes = []
    if os.path.exists(galigner.INDEXES_PATH):
        inside_indexes = os.listdir(galigner.INDEXES_PATH)
    if os.path.exists(os.path.join(MONICA_ROOT, '.outside_indexes')):
        outside_indexes = open(os.path.join(MONICA_ROOT, '.outside_indexes'), 'r').read().splitlines()

    # When pre-built indexes are provided
    if indexes:
        for index in indexes:
            if index in inside_indexes:
                index_folder = os.path.join(galigner.INDEXES_PATH, index)
            elif index in [os.path.basename(index) for index in outside_indexes]:
                index_folder = index
            else:
                index_folder = index
                with open(os.path.join(MONICA_ROOT, '.outside_indexes'), 'a+') as outside_indexes:
                    outside_indexes.seek(0)
                    data = outside_indexes.read(100)
                    if len(data) > 0:
                        outside_indexes.write("\n")
                    outside_indexes.write(index_folder)
            if os.path.basename(index_folder)[:2] == 'F_':
                focus_indexes_paths += [os.path.join(index_folder, file) for file in os.listdir(index_folder) if
                                        file.endswith('.mmi')]
            else:
                indexes_paths += [os.path.join(index_folder, file) for file in os.listdir(index_folder) if
                                  file.endswith('.mmi')]
            if os.path.basename(index)[0:2] == 'G_':
                index = os.path.basename(index)
                for i, j in zip(index.split(sep='_')[1:-2:2], index.split(sep='_')[2:-2:2]):
                    if j != 'm':
                        guest_species.append(i+'_'+j)
                    else:
                        guest_species.append(i)
            elif os.path.basename(index)[0:2] == 'H_':
                index = os.path.basename(index)
                host_species += [i + '_' + j for i, j in zip(index.split(sep='_')[1:-2:2], index.split(sep='_')[2:-2:2])]

    # Guests management: genomes retrieval, database and index formation
    if guests:
        print('Building guests databases and indexes...')
        guests_string = 'G_' + '_'.join(guests) + '_m_' + mode + '_im_' + str(
            round(indexing_memory * (10 ** -9), 3))
        guests = map(lambda guest: ' '.join(guest.split(sep='_')), guests)
        guests_databases_path = os.path.join(gdatabase.DATABASES_PATH, guests_string)

        # Genomes ftp selection and download:
        ftp_table = gfetcher.ftp_selector(mode=mode, species=guests)

        genomes, updated = gfetcher.fetcher(ftp_table, oldies_path=oldies_path, keep_genomes=keep_genomes,
                                            format_genomes=format_genomes)

        if guests_string not in inside_indexes and guests_string not in [os.path.basename(index) for index in
                                                                         outside_indexes] \
                and guests_string not in [os.path.basename(index) for index in indexes]:
            # Updated to see if genomes composing the database have new official versions
            # Database building and indexing:
            databases, genomes_length = gdatabase.multi_threaded_builder(genomes=genomes, max_chunk_size=max_chunk_size,
                                                                         databases_path=guests_databases_path,
                                                                         keep_genomes=keep_genomes, n_threads=n_threads)

            indexes_paths = indexes_paths + galigner.indexer(databases,
                                                             indexes_path=os.path.join(galigner.INDEXES_PATH,
                                                                                       guests_string))
        elif guests_string in [os.path.basename(index) for index in indexes]:
            index = [index for index in indexes if guests_string == os.path.basename(index)][-1]
            if updated:
                databases, genomes_length = gdatabase.multi_threaded_builder(genomes=genomes,
                                                                             max_chunk_size=max_chunk_size,
                                                                             databases_path=guests_databases_path,
                                                                             keep_genomes=keep_genomes,
                                                                             n_threads=n_threads)

                indexes_paths = indexes_paths + galigner.indexer(databases, indexes_path=index)
        elif guests_string in [os.path.basename(index) for index in outside_indexes]:
            index = [index for index in outside_indexes if guests_string == os.path.basename(index)][-1]
            if updated:
                databases, genomes_length = gdatabase.multi_threaded_builder(genomes=genomes,
                                                                             max_chunk_size=max_chunk_size,
                                                                             databases_path=guests_databases_path,
                                                                             keep_genomes=keep_genomes,
                                                                             n_threads=n_threads)

                indexes_paths = indexes_paths + galigner.indexer(databases, indexes_path=index)
            else:
                indexes_paths += [os.path.join(index, file) for file in os.listdir(index) if file.endswith('.mmi')]
        else:
            if updated:
                databases, genomes_length = gdatabase.multi_threaded_builder(genomes=genomes,
                                                                             max_chunk_size=max_chunk_size,
                                                                             databases_path=guests_databases_path,
                                                                             keep_genomes=keep_genomes,
                                                                             n_threads=n_threads)

                indexes_paths = indexes_paths + galigner.indexer(databases,
                                                                 indexes_path=os.path.join(galigner.INDEXES_PATH,
                                                                                           guests_string))
            else:
                index_folder = os.path.join(galigner.INDEXES_PATH, guests_string)
                indexes_paths += [os.path.join(index_folder, file) for file in os.listdir(index_folder) if
                                  file.endswith('.mmi')]

    # Hosts management: genomes retrieval, database and index formation
    if hosts:
        print('Building hosts databases and indexes...')
        hosts = map(lambda host: ' '.join(host.split(sep='_')), hosts)
        hosts_ftp_table = gfetcher.ftp_selector(mode='single', species=hosts)
        #    host_ftp_table = ftp_table.append(host_ftp_table, ignore_index=True)
        host_genomes, host_updated = gfetcher.fetcher(hosts_ftp_table, oldies_path=oldies_path,
                                                      keep_genomes=keep_genomes,
                                                      format_genomes=format_genomes)
        for genome in host_genomes:
            host_string = 'H_' + genome[1][0] + '_im_' + str(round(indexing_memory * (10 ** -9), 3))
            host_databases_path = os.path.join(gdatabase.DATABASES_PATH, host_string)
            if host_string not in inside_indexes and host_string not in [os.path.basename(index) for index in
                                                                         outside_indexes] \
                    and host_string not in [os.path.basename(index) for index in indexes]:
                databases, genomes_length = gdatabase.multi_threaded_builder(genomes=[genome],
                                                                             max_chunk_size=max_chunk_size,
                                                                             databases_path=host_databases_path,
                                                                             keep_genomes=keep_genomes,
                                                                             n_threads=n_threads)
                indexes_paths = indexes_paths + galigner.indexer(databases,
                                                                 indexes_path=os.path.join(galigner.INDEXES_PATH,
                                                                                           host_string))
            elif host_string in [os.path.basename(index) for index in indexes]:
                index = [index for index in indexes if host_string == os.path.basename(index)][-1]
                if host_updated:
                    databases, genomes_length = gdatabase.multi_threaded_builder(genomes=[genome],
                                                                                 max_chunk_size=max_chunk_size,
                                                                                 databases_path=host_databases_path,
                                                                                 keep_genomes=keep_genomes,
                                                                                 n_threads=n_threads)

                    indexes_paths = indexes_paths + galigner.indexer(databases, indexes_path=index)
            elif host_string in [os.path.basename(index) for index in outside_indexes]:
                index = [index for index in outside_indexes if host_string == os.path.basename(index)][-1]
                if host_updated:
                    databases, genomes_length = gdatabase.multi_threaded_builder(genomes=[genome],
                                                                                 max_chunk_size=max_chunk_size,
                                                                                 databases_path=host_databases_path,
                                                                                 keep_genomes=keep_genomes,
                                                                                 n_threads=n_threads)

                    indexes_paths = indexes_paths + galigner.indexer(databases, indexes_path=index)
                else:
                    indexes_paths += [os.path.join(index, file) for file in os.listdir(index) if file.endswith('.mmi')]
            else:
                if host_updated:
                    databases, genomes_length = gdatabase.multi_threaded_builder(genomes=[genome],
                                                                                 max_chunk_size=max_chunk_size,
                                                                                 databases_path=host_databases_path,
                                                                                 keep_genomes=keep_genomes,
                                                                                 n_threads=n_threads)

                    indexes_paths = indexes_paths + galigner.indexer(databases,
                                                                     indexes_path=os.path.join(galigner.INDEXES_PATH,
                                                                                               host_string))
                else:
                    index_folder = os.path.join(galigner.INDEXES_PATH, host_string)
                    indexes_paths += [os.path.join(index_folder, file) for file in os.listdir(index_folder) if
                                      file.endswith('.mmi')]

    # Focus mode management: genomes retrieval, database and index formation
    if focus_species:
        print('Building focus mode databases and indexes...')
        focus_string = 'F_' + '_'.join(focus_species) + '_im_' + str(round(indexing_memory * (10 ** -9), 3))
        focus_species = list(map(lambda specie: ' '.join(specie.split(sep='_')), focus_species))
        focus_databases_path = os.path.join(gdatabase.DATABASES_PATH, focus_string)
        focus_ftp_table = gfetcher.ftp_selector(mode='focus', species=focus_species)
        focus_genomes, focus_updated = gfetcher.focus_fetcher(focus_ftp_table, oldies_path=oldies_path,
                                                              keep_genomes=keep_genomes)

        # Focus mode management: genomes retrieval, database and index formation
        if focus_string not in inside_indexes and focus_string not in [os.path.basename(index) for index in
                                                                       outside_indexes] \
                and focus_string not in [os.path.basename(index) for index in indexes]:
            focus_databases, focus_genomes_length = gdatabase.multi_threaded_builder(genomes=focus_genomes,
                                                                                     max_chunk_size=max_chunk_size,
                                                                                     keep_genomes=keep_genomes,
                                                                                     n_threads=n_threads,
                                                                                     databases_path=focus_databases_path)
            focus_indexes_paths += galigner.indexer(focus_databases,
                                                    indexes_path=os.path.join(galigner.INDEXES_PATH, focus_string))
        elif focus_string in [os.path.basename(index) for index in indexes]:
            index = [index for index in indexes if focus_string == os.path.basename(index)][-1]
            if focus_updated:
                focus_databases, focus_genomes_length = gdatabase.multi_threaded_builder(genomes=focus_genomes,
                                                                                         max_chunk_size=max_chunk_size,
                                                                                         keep_genomes=keep_genomes,
                                                                                         n_threads=n_threads,
                                                                                         databases_path=focus_databases_path)
                focus_indexes_paths += galigner.indexer(focus_databases, indexes_path=index)
        elif focus_string in [os.path.basename(index) for index in outside_indexes]:
            index = [index for index in outside_indexes if focus_string == os.path.basename(index)][-1]
            if focus_updated:
                focus_databases, focus_genomes_length = gdatabase.multi_threaded_builder(genomes=focus_genomes,
                                                                                         max_chunk_size=max_chunk_size,
                                                                                         keep_genomes=keep_genomes,
                                                                                         n_threads=n_threads,
                                                                                         databases_path=focus_databases_path)

                focus_indexes_paths += galigner.indexer(focus_databases, indexes_path=index)
            else:
                focus_indexes_paths += [os.path.join(index, file) for file in os.listdir(index) if
                                        file.endswith('.mmi')]
        else:
            if focus_updated:
                focus_databases, focus_genomes_length = gdatabase.multi_threaded_builder(genomes=focus_genomes,
                                                                                         max_chunk_size=max_chunk_size,
                                                                                         keep_genomes=keep_genomes,
                                                                                         n_threads=n_threads,
                                                                                         databases_path=focus_databases_path)
                focus_indexes_paths += galigner.indexer(focus_databases,
                                                        indexes_path=os.path.join(galigner.INDEXES_PATH, focus_string))
            else:
                index_folder = os.path.join(galigner.INDEXES_PATH, focus_string)
                focus_indexes_paths += [os.path.join(index_folder, file) for file in os.listdir(index_folder) if
                                        file.endswith('.mmi')]

    with open(os.path.join(gfetcher.GENOMES_PATH, 'going_to_enter_alignment'), 'wb'):
        pass

    if indexes_paths:
        if indexes and focus_indexes_paths:
            focus_species = list(map(lambda specie: '_'.join(specie.split(sep=' ')), focus_species))
            for index in indexes:
                index = os.path.basename(index)
                if index[:2] == 'F_':
                    focus_species += [i + '_' + j for i, j in
                                      zip(index.split(sep='_')[1:-2:2], index.split(sep='_')[2:-2:2])]

        # Alignment and normalization:
        alignment = galigner.multi_threaded_aligner(input_folder, indexes_paths, mode=alignment_mode,
                                                    n_threads=n_threads,
                                                    focus_species=focus_species, output_folder=output_folder)

        if alignment:

            raw_alignment_df = galigner.alignment_to_data_frame(alignment, output_folder=output_folder,
                                                                filename='raw_monica.dataframe')

            norm_alignment = galigner.normalizer(alignment)

            norm_alignment_df = galigner.alignment_to_data_frame(norm_alignment, output_folder=output_folder)

            # Plotting:
            barplot.plotter(norm_alignment_df, raw_alignment_df, output_folder=output_folder, palette='jet',
                            reads_threshold=reads_threshold, hosts=host_species, guests=guest_species, mode=mode,
                            show_legend=show_legend, auto_open=auto_open_plot)

            if focus_indexes_paths:
                print('Entering focus mode')
                focus_input_folder = os.path.join(input_folder, 'focus')
                focus_output_folder = os.path.join(output_folder, 'focus')
                if not os.path.exists(focus_output_folder):
                    os.mkdir(focus_output_folder)

                helpers.initializer(focus_output_folder)

                focus_alignment = galigner.multi_threaded_aligner(focus_input_folder, focus_indexes_paths,
                                                                  mode=alignment_mode,
                                                                  n_threads=n_threads,
                                                                  output_folder=focus_output_folder)
                if focus_alignment:
                    focus_raw_alignment_df = galigner.alignment_to_data_frame(focus_alignment,
                                                                              output_folder=focus_output_folder,
                                                                              filename='raw_monica.dataframe')
                    focus_norm_alignment = galigner.normalizer(focus_alignment)
                    focus_norm_alignment_df = galigner.alignment_to_data_frame(focus_norm_alignment,
                                                                               output_folder=focus_output_folder)
                    barplot.plotter(focus_norm_alignment_df, focus_raw_alignment_df, output_folder=focus_output_folder,
                                    palette='jet', reads_threshold=0, guests=focus_species, mode='focus',
                                    show_legend=show_legend, auto_open=auto_open_plot)


def main_build_index(args):
    # Input handling:
    with open(os.path.join(MONICA_ROOT, 'monica.params'), 'w') as params:
        params.write(str(args))

    n_threads = args.threads

    if not args.focus_species:
        focus_species = []
    else:
        focus_species = args.focus_species

    if args.indexing_memory:
        indexing_memory = args.indexing_memory
    else:
        indexing_memory = psutil.virtual_memory().total / 4

    max_chunk_size = indexing_memory / 16

    mode = args.mode
    hosts = args.host_specie
    guests = args.guest_species

    if args.indexes:
        indexes_folder = args.indexes[0]
    else:
        indexes_folder = galigner.INDEXES_PATH

    if args.keep_genomes == 'yes':
        keep_genomes = True
    else:
        keep_genomes = False
    oldies_path = args.genomes_folder
    format_genomes = args.format_genomes

    indexes_paths = []
    focus_indexes_paths = []
    inside_indexes = []
    outside_indexes = []
    if os.path.exists(galigner.INDEXES_PATH):
        inside_indexes = os.listdir(galigner.INDEXES_PATH)
    if os.path.exists(os.path.join(MONICA_ROOT, '.outside_indexes')):
        outside_indexes = open(os.path.join(MONICA_ROOT, '.outside_indexes'), 'r').read().splitlines()

    # Guests management: genomes retrieval, database and index formation
    if guests:
        print('Building guests databases and indexes...')
        guests_string = 'G_' + '_'.join(guests) + '_m_' + mode + '_im_' + str(round(indexing_memory * (10 ** -9), 3))
        guests = map(lambda guest: ' '.join(guest.split(sep='_')), guests)
        guests_databases_path = os.path.join(gdatabase.DATABASES_PATH, guests_string)

        # Genomes ftp selection and download:
        ftp_table = gfetcher.ftp_selector(mode=mode, species=guests)

        genomes, updated = gfetcher.fetcher(ftp_table, oldies_path=oldies_path, keep_genomes=keep_genomes,
                                            format_genomes=format_genomes)

        if guests_string not in inside_indexes and guests_string not in [os.path.basename(index) for index in
                                                                         outside_indexes]:
            # Updated to see if genomes composing the database have new official versions
            # Database building and indexing:
            databases, genomes_length = gdatabase.multi_threaded_builder(genomes=genomes, max_chunk_size=max_chunk_size,
                                                                         databases_path=guests_databases_path,
                                                                         keep_genomes=keep_genomes, n_threads=n_threads)

            indexes_paths = indexes_paths + galigner.indexer(databases,
                                                             indexes_path=os.path.join(indexes_folder, guests_string))

            if indexes_folder != galigner.INDEXES_PATH:
                with open(os.path.join(MONICA_ROOT, '.outside_indexes'), 'a+') as outside_indexes:
                    outside_indexes.seek(0)
                    data = outside_indexes.read(100)
                    if len(data) > 0:
                        outside_indexes.write("\n")
                    outside_indexes.write(os.path.join(indexes_folder, guests_string))

        elif guests_string in [os.path.basename(index) for index in outside_indexes]:
            index = [index for index in outside_indexes if guests_string == os.path.basename(index)][-1]
            if updated:
                databases, genomes_length = gdatabase.multi_threaded_builder(genomes=genomes,
                                                                             max_chunk_size=max_chunk_size,
                                                                             databases_path=guests_databases_path,
                                                                             keep_genomes=keep_genomes,
                                                                             n_threads=n_threads)

                indexes_paths = indexes_paths + galigner.indexer(databases, indexes_path=index)
            else:
                indexes_paths += [os.path.join(index, file) for file in os.listdir(index) if file.endswith('.mmi')]
            if indexes_folder != index:
                print('Selected index is already present on the machine. Cloning it to specified location')
                shutil.copytree(index, os.path.join(indexes_folder, guests_string))

        else:
            if updated:
                databases, genomes_length = gdatabase.multi_threaded_builder(genomes=genomes,
                                                                             max_chunk_size=max_chunk_size,
                                                                             databases_path=guests_databases_path,
                                                                             keep_genomes=keep_genomes,
                                                                             n_threads=n_threads)

                indexes_paths = indexes_paths + galigner.indexer(databases,
                                                                 indexes_path=os.path.join(galigner.INDEXES_PATH,
                                                                                           guests_string))
            else:
                index_folder = os.path.join(galigner.INDEXES_PATH, guests_string)
                indexes_paths += [os.path.join(index_folder, file) for file in os.listdir(index_folder) if
                                  file.endswith('.mmi')]
            if indexes_folder != galigner.INDEXES_PATH:
                print('Selected index is already present on the machine. Cloning it to specified location')
                index = os.path.join(galigner.INDEXES_PATH, guests_string)
                shutil.copytree(index, os.path.join(indexes_folder, guests_string))

    # Hosts management: genomes retrieval, database and index formation
    if hosts:
        print('Building hosts databases and indexes...')
        hosts = map(lambda host: ' '.join(host.split(sep='_')), hosts)
        hosts_ftp_table = gfetcher.ftp_selector(mode='single', species=hosts)
        #    host_ftp_table = ftp_table.append(host_ftp_table, ignore_index=True)
        host_genomes, host_updated = gfetcher.fetcher(hosts_ftp_table, oldies_path=oldies_path,
                                                      keep_genomes=keep_genomes,
                                                      format_genomes=format_genomes)
        for genome in host_genomes:
            host_string = 'H_' + genome[1][0] + '_im_' + str(round(indexing_memory * (10 ** -9), 3))
            host_databases_path = os.path.join(gdatabase.DATABASES_PATH, host_string)
            if host_string not in inside_indexes and host_string not in [os.path.basename(index) for index in
                                                                         outside_indexes]:
                # Updated to see if genomes composing the database have new official versions
                # Database building and indexing:
                databases, genomes_length = gdatabase.multi_threaded_builder(genomes=[genome],
                                                                             max_chunk_size=max_chunk_size,
                                                                             databases_path=host_databases_path,
                                                                             keep_genomes=keep_genomes,
                                                                             n_threads=n_threads)

                indexes_paths = indexes_paths + galigner.indexer(databases,
                                                                 indexes_path=os.path.join(indexes_folder,
                                                                                           host_string))

                if indexes_folder != galigner.INDEXES_PATH:
                    with open(os.path.join(MONICA_ROOT, '.outside_indexes'), 'a+') as outside_indexes:
                        outside_indexes.seek(0)
                        data = outside_indexes.read(100)
                        if len(data) > 0:
                            outside_indexes.write("\n")
                        outside_indexes.write(os.path.join(indexes_folder, host_string))

            elif host_string in [os.path.basename(index) for index in outside_indexes]:
                index = [index for index in outside_indexes if host_string == os.path.basename(index)][-1]
                if host_updated:
                    databases, genomes_length = gdatabase.multi_threaded_builder(genomes=[genome],
                                                                                 max_chunk_size=max_chunk_size,
                                                                                 databases_path=host_databases_path,
                                                                                 keep_genomes=keep_genomes,
                                                                                 n_threads=n_threads)

                    indexes_paths = indexes_paths + galigner.indexer(databases, indexes_path=index)
                else:
                    indexes_paths += [os.path.join(index, file) for file in os.listdir(index) if file.endswith('.mmi')]
                if indexes_folder != index:
                    print('Selected index is already present on the machine. Cloning it to specified location')
                    shutil.copytree(index, os.path.join(indexes_folder, host_string))

            else:
                if host_updated:
                    databases, genomes_length = gdatabase.multi_threaded_builder(genomes=[genome],
                                                                                 max_chunk_size=max_chunk_size,
                                                                                 databases_path=host_databases_path,
                                                                                 keep_genomes=keep_genomes,
                                                                                 n_threads=n_threads)

                    indexes_paths = indexes_paths + galigner.indexer(databases,
                                                                     indexes_path=os.path.join(galigner.INDEXES_PATH,
                                                                                               host_string))
                else:
                    index_folder = os.path.join(galigner.INDEXES_PATH, host_string)
                    indexes_paths += [os.path.join(index_folder, file) for file in os.listdir(index_folder) if
                                      file.endswith('.mmi')]
                if indexes_folder != galigner.INDEXES_PATH:
                    print('Selected index is already present on the machine. Cloning it to specified location')
                    index = os.path.join(galigner.INDEXES_PATH, host_string)
                    shutil.copytree(index, os.path.join(indexes_folder, host_string))

    # Focus mode management: genomes retrieval, database and index formation
    if focus_species:
        print('Building focus mode databases and indexes...')
        focus_string = 'F_' + '_'.join(focus_species) + '_im_' + str(round(indexing_memory * (10 ** -9), 3))
        focus_databases_path = os.path.join(gdatabase.DATABASES_PATH, focus_string)
        focus_species = map(lambda specie: ' '.join(specie.split(sep='_')), focus_species)
        focus_ftp_table = gfetcher.ftp_selector(mode='focus', species=focus_species)
        focus_genomes, focus_updated = gfetcher.focus_fetcher(focus_ftp_table, oldies_path=oldies_path,
                                                              keep_genomes=keep_genomes)

        if focus_string not in inside_indexes and focus_string not in [os.path.basename(index) for index in
                                                                       outside_indexes]:
            focus_databases, focus_genomes_length = gdatabase.multi_threaded_builder(genomes=focus_genomes,
                                                                                     max_chunk_size=max_chunk_size,
                                                                                     keep_genomes=keep_genomes,
                                                                                     n_threads=n_threads,
                                                                                     databases_path=focus_databases_path)

            focus_indexes_paths += galigner.indexer(focus_databases,
                                                    indexes_path=os.path.join(indexes_folder, focus_string))

            if indexes_folder != galigner.INDEXES_PATH:
                with open(os.path.join(MONICA_ROOT, '.outside_indexes'), 'a+') as outside_indexes:
                    outside_indexes.seek(0)
                    data = outside_indexes.read(100)
                    if len(data) > 0:
                        outside_indexes.write("\n")
                    outside_indexes.write(os.path.join(indexes_folder, focus_string))

        elif focus_string in [os.path.basename(index) for index in outside_indexes]:
            index = [index for index in outside_indexes if focus_string == os.path.basename(index)][-1]
            if focus_updated:
                focus_databases, focus_genomes_length = gdatabase.multi_threaded_builder(genomes=focus_genomes,
                                                                                         max_chunk_size=max_chunk_size,
                                                                                         databases_path=focus_databases_path,
                                                                                         keep_genomes=keep_genomes,
                                                                                         n_threads=n_threads)

                focus_indexes_paths += galigner.indexer(focus_databases, indexes_path=index)
            else:
                focus_indexes_paths += [os.path.join(index, file) for file in os.listdir(index) if
                                        file.endswith('.mmi')]
            if indexes_folder != index:
                print('Selected index is already present on the machine. Cloning it to specified location')
                shutil.copytree(index, os.path.join(indexes_folder, focus_string))

        else:
            if focus_updated:
                focus_databases, focus_genomes_length = gdatabase.multi_threaded_builder(genomes=focus_genomes,
                                                                                         max_chunk_size=max_chunk_size,
                                                                                         databases_path=focus_databases_path,
                                                                                         keep_genomes=keep_genomes,
                                                                                         n_threads=n_threads)

                focus_indexes_paths += galigner.indexer(focus_databases,
                                                        indexes_path=os.path.join(galigner.INDEXES_PATH, focus_string))
            else:
                index_folder = os.path.join(galigner.INDEXES_PATH, focus_string)
                focus_indexes_paths += [os.path.join(index_folder, file) for file in os.listdir(index_folder) if
                                        file.endswith('.mmi')]
            if indexes_folder != galigner.INDEXES_PATH:
                print('Selected index is already present on the machine. Cloning it to specified location')
                index = os.path.join(galigner.INDEXES_PATH, focus_string)
                shutil.copytree(index, os.path.join(indexes_folder, focus_string))
    print('\nBuilt the following indexes:\n{}\n{}'.format('\n'.join(indexes_paths), '\n'.join(focus_indexes_paths)))


def main_list_indexes(args):
    indexes = []
    outside_indexes = []
    indexes += [folder for folder in os.listdir(galigner.INDEXES_PATH) if
                os.path.isdir(os.path.join(galigner.INDEXES_PATH, folder))]
    if os.path.exists(os.path.join(MONICA_ROOT, '.outside_indexes')):
        outside_indexes = open(os.path.join(MONICA_ROOT, '.outside_indexes'), 'r').read().splitlines()
    indexes += [os.path.basename(folder) for folder in outside_indexes]
    guests_indexes = [index for index in indexes if index[:2] == 'G_']
    hosts_indexes = [index for index in indexes if index[:2] == 'H_']
    focus_indexes = [index for index in indexes if index[:2] == 'F_']
    print('\nGuest indexes available:')
    print('\t' + '\n\t'.join(guests_indexes))
    print('\n\nHosts indexes available:')
    print('\t' + '\n\t'.join(hosts_indexes))
    print('\n\nFocus indexes available:')
    print('\t' + '\n\t'.join(focus_indexes))


def main_from_plotting(args):
    # Input handling
    data_frame_path = args.data_frame
    ex_out_dir = os.path.dirname(data_frame_path)
    r_data_frame_path = os.path.join(ex_out_dir, 'raw_monica.dataframe')

    if args.output_folder:
        output_folder = args.output_folder
    else:
        output_folder = os.path.join(ex_out_dir, 'monica_output')
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    auto_open_plot = not args.not_auto_open_plot
    show_legend = not args.not_show_legend
    reads_threshold = args.reads_threshold

    if 'monica.params' in os.listdir(ex_out_dir):
        with open(os.path.join(ex_out_dir, 'monica.params'), 'r') as ex_params:
            pass
    else:
        # Importing data frame
        n_alignment_df = pd.read_csv(data_frame_path, index_col=(0, 1))
        r_alignment_df = pd.read_csv(r_data_frame_path, index_col=(0, 1))
        # Plotting
        barplot.plotter(n_alignment_df, r_alignment_df, output_folder=output_folder, palette='jet',
                        reads_threshold=reads_threshold, show_legend=show_legend, auto_open=auto_open_plot)


def main_initialize(args):
    # To initialize monica in a given root
    if args.monica_root_folder != MONICA_ROOT:
        monica_root = os.path.join(args.monica_root_folder, '.monica')
        with open(os.path.join(MONICA_ROOT, '.root'), 'w') as root:
            root.write(monica_root)
        if not os.path.exists(monica_root):
            os.mkdir(monica_root)
        print(f'monica root moved to {monica_root}')


if __name__ == '__main__':
    main()
