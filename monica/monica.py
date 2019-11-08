#! /usr/bin/env python3.6

# TODO log system when launching monica from scratch to save all parameters, to keep host, guests and mode

import os
import argparse

import psutil
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
    parser.add_argument('-t', '--threads', type=int, default=3,
                        help='The number of threads to be used')
    parser.add_argument('-a', '--alignment_mode', default='query_length',
                              help='The alignment mode to be performed')
    parser.add_argument('-g', '--genomes_folder', default=gfetcher.OLDIES_PATH,
                              help='The folder where already stored genomes will be found '
                                   'and new genomes will be stored if --keep_genomes option is set to "yes". '
                                   'Default is {}'.format(gfetcher.OLDIES_PATH))
    parser.add_argument('-k', '--keep_genomes', choices=['yes', 'no'], default='yes',
                              help='Choose if genomes are to be kept in -g. Default is "yes"')
    parser.add_argument('-im', '--indexing_memory', type=helpers.human_readable,
                              help='The maximum memory the user wishes to be used during indexing. '
                                   'Must end with a unit of measure among B|K|M|G|T.'
                                   'Default is half of total memory')
    parser.add_argument('--not_auto_open_plot', action='store_true',
                        help='If wishing NOT to plot the results at the and of the analysis')
    parser.add_argument('--not_show_legend', action='store_true',
                        help='If wishing NOT to show legend')
    parser.add_argument('-r', '--reads_threshold', default=0,
                        help='Threshold of reads for a specie in each sample that will be displayed in the plot. '
                             'If lower in all samples species not displayed. '
                             'Default is {}'.format(barplot.READS_THRESHOLD))

    subparsers = parser.add_subparsers()

    from_scratch = subparsers.add_parser('from_scratch', help='To run monica from the beginning')
    from_scratch.add_argument('-H', '--host_specie',
                        help='Host species where samples come from, underscore separated ("Genus_specie")')
    from_scratch.add_argument('-G', '--guest_species', nargs='*',
                        help='Guest species hypothesized to be present in the sample, '
                             'underscore separated ("Genus_specie"), but also higher levels of the tree of life, '
                             'like order or genus only')
    from_scratch.add_argument('-F', '--focus_species', nargs='*',
                        help='Species of interest among the whole guests selected, focus analysis on their subspecies')
    from_scratch.add_argument('--format_genomes',
                        help='Folder, call it together with -m and -H/-G, when assembly-formatted genomes '
                             'are present in it without being formatted for monica')
    from_scratch.add_argument('-m', '--mode', choices=['single', 'all', 'overnight'],
                              help="The mode of monica's execution: single, all or overnight")
    from_scratch.set_defaults(func=main_from_scratch)

    from_alignment = subparsers.add_parser('from_alignment', help='To run monica given an already built index/database')
    from_alignment.add_argument('-F', '--focus_species', nargs='*',
                        help='Species of interest among the whole guests selected, focus analysis on their subspecies')
    from_alignment.add_argument('-upi', '--use_prebuilt_indexes', action='store_true',
                        help='Use the index of the previous analysis, useful for frequent use of the same genomes')
    from_alignment.add_argument('-upd', '--use_prebuilt_databases', action='store_true')
    from_alignment.add_argument('-H', '--host_specie',
                        help='Host species where samples come from, underscore separated ("Genus_specie")')
    from_alignment.add_argument('-G', '--guest_species', nargs='*',
                        help='Guest species hypothesized to be present in the sample, '
                             'underscore separated ("Genus_specie"), but also higher levels of the tree of life, '
                             'like order or genus only')
    from_alignment.add_argument('-m', '--mode', choices=['single', 'all', 'overnight'],
                        help="The mode of monica's execution: single, all or overnight")
    from_alignment.set_defaults(func=main_from_alignment)

    plot_only = subparsers.add_parser('plot_only', help='To generate the plot only, given a table of monica results')
    plot_only.add_argument('-d', '--data_frame', help='The monica result data frame to generate the plot from')
    plot_only.set_defaults(func=main_from_plotting)

    args = parser.parse_args()
    args.func(args)
    return args


def main_from_scratch(args):

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
    with open(os.path.join(output_folder, 'monica.params'), 'w') as params:
        params.write(str(args))

    n_threads = args.threads
    alignment_mode = args.alignment_mode

    focus_species = args.focus_species

    auto_open_plot = not args.not_auto_open_plot
    show_legend = not args.not_show_legend
    reads_threshold = args.reads_threshold

    if args.indexing_memory:
        indexing_memory = args.indexing_memory
    else:
        indexing_memory = psutil.virtual_memory().total/4

    max_chunk_size = indexing_memory/16

    mode = args.mode
    host = args.host_specie
    guests = args.guest_species

    helpers.initializer(output_folder)

    if guests:
        guests = map(lambda guest: ' '.join(guest.split(sep='_')), guests)
    if args.keep_genomes == 'yes':
        keep_genomes = True
    else:
        keep_genomes = False
    oldies_path = args.genomes_folder
    format_genomes = args.format_genomes

    # Genomes ftp selection and download:
    ftp_table = gfetcher.ftp_selector(mode=mode, species=guests)

    if host:
        host = ' '.join(host.split(sep='_'))
        host_ftp_table = gfetcher.ftp_selector(mode='single', species=[host])
        ftp_table = ftp_table.append(host_ftp_table, ignore_index=True)

    genomes = gfetcher.fetcher(ftp_table, oldies_path=oldies_path, keep_genomes=keep_genomes,
                               format_genomes=format_genomes)

    # Database building and indexing:
    databases, genomes_length = gdatabase.multi_threaded_builder(genomes=genomes, max_chunk_size=max_chunk_size,
                                                                 databases_path=gdatabase.DATABASES_PATH,
                                                                 keep_genomes=keep_genomes, n_threads=n_threads)

    indexes_paths = galigner.indexer(databases)

    # Focus mode indexing building block:
    if focus_species:
        focus_input_folder = os.path.join(input_folder, 'focus')
        focus_output_folder = os.path.join(output_folder, 'focus')
        if not os.path.exists(focus_output_folder):
            os.mkdir(focus_output_folder)

        helpers.initializer(focus_output_folder)

        focus_databases_path = os.path.join(gdatabase.DATABASES_PATH, 'focus')
        focus_species = map(lambda specie: ' '.join(specie.split(sep='_')), focus_species)
        focus_ftp_table = gfetcher.ftp_selector(mode='focus', species=focus_species)
        focus_genomes = gfetcher.focus_fetcher(focus_ftp_table, oldies_path=oldies_path, keep_genomes=keep_genomes)
        focus_databases, focus_genomes_length = gdatabase.multi_threaded_builder(genomes=focus_genomes,
                                                                                 max_chunk_size=max_chunk_size,
                                                                                 keep_genomes=keep_genomes,
                                                                                 n_threads=n_threads,
                                                                                 databases_path=focus_databases_path)
        focus_indexes_paths = galigner.indexer(focus_databases, indexes_path=os.path.join(galigner.INDEXES_PATH, 'focus'))

    with open(os.path.join(gfetcher.GENOMES_PATH, 'going_to_enter_alignment'), 'wb'):
        pass

    # Alignment and normalization:
    alignment = galigner.multi_threaded_aligner(input_folder, indexes_paths, mode=alignment_mode, n_threads=n_threads,
                                                focus_species=args.focus_species, output_folder=output_folder)

    raw_alignment_df = galigner.alignment_to_data_frame(alignment, output_folder=output_folder,
                                                        filename='raw_monica.dataframe')

    norm_alignment = galigner.normalizer(alignment, genomes_length)

    norm_alignment_df = galigner.alignment_to_data_frame(norm_alignment, output_folder=output_folder)

    # Plotting:
    barplot.plotter(norm_alignment_df, raw_alignment_df, output_folder=output_folder, palette='jet',
                    reads_threshold=reads_threshold, host=host, guests=args.guest_species, mode=mode,
                    show_legend=show_legend, auto_open=auto_open_plot)

    if focus_species:
        focus_alignment = galigner.multi_threaded_aligner(focus_input_folder, focus_indexes_paths, mode=alignment_mode,
                                                          n_threads=n_threads, output_folder=focus_output_folder)
        focus_raw_alignment_df = galigner.alignment_to_data_frame(focus_alignment, output_folder=focus_output_folder,
                                                                  filename='raw_monica.dataframe')
        focus_norm_alignment = galigner.normalizer(focus_alignment, focus_genomes_length)
        focus_norm_alignment_df = galigner.alignment_to_data_frame(focus_norm_alignment, 
                                                                   output_folder=focus_output_folder)
        barplot.plotter(focus_norm_alignment_df, focus_raw_alignment_df, output_folder=focus_output_folder, 
                        palette='jet', reads_threshold=0, guests=args.focus_species, mode='focus',
                        show_legend=show_legend, auto_open=auto_open_plot)
        
        
def main_from_alignment(args):

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
    with open(os.path.join(output_folder, 'monica.params'), 'w') as params:
        params.write(str(args))

    focus_species = args.focus_species

    n_threads = args.threads
    alignment_mode = args.alignment_mode

    auto_open_plot = not args.not_auto_open_plot
    show_legend = not args.not_show_legend
    reads_threshold = args.reads_threshold

    mode = args.mode
    host = args.host_specie
    guests = args.guest_species

    helpers.initializer(output_folder)

    use_prebuilt_databases = args.use_prebuilt_databases
    use_prebuilt_indexes = args.use_prebuilt_indexes

    if use_prebuilt_databases:
        databases = os.path.join(gdatabase.DATABASES_PATH)

        indexes = galigner.indexer(databases)

    elif use_prebuilt_indexes:
        indexes_paths = [file for file in os.listdir(galigner.INDEXES_PATH)]

    with open(os.path.join(gfetcher.GENOMES_PATH, 'going_to_enter_alignment'), 'wb'):
        pass

    # Alignment and normalization:
    alignment = galigner.multi_threaded_aligner(input_folder, indexes_paths, mode=alignment_mode, n_threads=n_threads,
                                                focus_species=args.focus_species, output_folder=output_folder)

    raw_alignment_df = galigner.alignment_to_data_frame(alignment, output_folder=output_folder,
                                                        filename='raw_monica.dataframe')

    norm_alignment = galigner.normalizer(alignment, databases_path=gdatabase.DATABASES_PATH)

    norm_alignment_df = galigner.alignment_to_data_frame(norm_alignment, output_folder=output_folder)

    # Plotting:
    barplot.plotter(norm_alignment_df, raw_alignment_df, output_folder=output_folder, palette='jet',
                    reads_threshold=reads_threshold, host=host, guests=guests, mode=mode,
                    show_legend=show_legend, auto_open=auto_open_plot)

    if focus_species:

        oldies_path = args.genomes_folder
        if args.keep_genomes == 'yes':
            keep_genomes = True
        else:
            keep_genomes = False
        if args.indexing_memory:
            indexing_memory = args.indexing_memory
        else:
            indexing_memory = psutil.virtual_memory().total / 4

        max_chunk_size = indexing_memory / 16
        focus_input_folder = os.path.join(input_folder, 'focus')
        focus_output_folder = os.path.join(output_folder, 'focus')
        if not os.path.exists(focus_output_folder):
            os.mkdir(focus_output_folder)

        helpers.initializer(focus_output_folder)

        focus_databases_path = os.path.join(gdatabase.DATABASES_PATH, 'focus')
        focus_species = map(lambda specie: ' '.join(specie.split(sep='_')), focus_species)
        focus_ftp_table = gfetcher.ftp_selector(mode='focus', species=focus_species)
        focus_genomes = gfetcher.focus_fetcher(focus_ftp_table, oldies_path=oldies_path, keep_genomes=keep_genomes)
        focus_databases, focus_genomes_length = gdatabase.multi_threaded_builder(genomes=focus_genomes,
                                                                                 max_chunk_size=max_chunk_size,
                                                                                 keep_genomes=keep_genomes,
                                                                                 n_threads=n_threads,
                                                                                 databases_path=focus_databases_path)
        focus_indexes = galigner.indexer(focus_databases, indexes_path=os.path.join(galigner.INDEXES_PATH, 'focus'))
        focus_alignment = galigner.multi_threaded_aligner(focus_input_folder, focus_indexes, mode=alignment_mode,
                                                          n_threads=n_threads, output_folder=focus_output_folder)
        focus_raw_alignment_df = galigner.alignment_to_data_frame(focus_alignment, output_folder=focus_output_folder,
                                                                  filename='raw_monica.dataframe')
        focus_norm_alignment = galigner.normalizer(focus_alignment, focus_genomes_length)
        focus_norm_alignment_df = galigner.alignment_to_data_frame(focus_norm_alignment,
                                                                   output_folder=focus_output_folder)
        barplot.plotter(focus_norm_alignment_df, focus_raw_alignment_df, output_folder=focus_output_folder,
                        palette='jet', reads_threshold=0, guests=args.focus_species, mode='focus',
                        show_legend=show_legend, auto_open=auto_open_plot)


def main_indexes_building_only():
    pass


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


if __name__ == '__main__':

    main()


