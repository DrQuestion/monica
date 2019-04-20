#! /home/drq/.venvs/monica/bin/python

import os
import argparse

import genomes.fetcher as gfetcher
import genomes.database as gdatabase
import genomes.aligner as galigner
import plots.barplot as barplot


def check_input():
    parser = argparse.ArgumentParser(prog='monica')

    input_folder = parser.add_mutually_exclusive_group(required=True)
    input_folder.add_argument('-q', '--query_folder',
                              help='The folder where input fastq are stored')
    input_folder.add_argument('-f5', '--fast5_folder',
                              help='The folder where input fast5 are stored (real-time use only)')
    parser.add_argument('-o', '--output_folder',
                        help='The folder where output files will be stored')
    parser.add_argument('-H', '--host_specie',
                        help='Host species where samples come from, underscore separated ("Genus_specie")')
    parser.add_argument('-G', '--guest_species', nargs='*',
                        help='Guest species hypothesized to be present in sample, underscore separated ("Genus_specie")')
    parser.add_argument('-k', '--keep_genomes', choices=['yes', 'no'], required=True,
                        help='Choose if genomes are to be kept in -g')
    parser.add_argument('-g', '--genomes_folder', default=gfetcher.OLDIES_PATH,
                        help='The folder where already stored genomes will be found '
                             'and new genomes will be stored if --keep_genomes option is set to "yes"')
    parser.add_argument('-t', '--threads', type=int, default=3,
                        help='The number of threads to be used')
    parser.add_argument('-m', '--mode', choices=['single', 'all', 'overnight'],
                        help="The mode of monica's execution: single, all or overnight")
    parser.add_argument('--format_genomes',
                        help='Folder, call it together with -m and -H/-G, when assembly-formatted genomes'
                             'are present in it without being formatted for monica')

    args = parser.parse_args()
    return args


def main():

    # Input handling
    args = check_input()
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
    host = args.host_specie
    guests = args.guest_species
    if guests:
        guests=map(lambda guest: ' '.join(guest.split(sep='_')), guests)
    if args.keep_genomes == 'yes':
        keep_genomes = True
    else:
        keep_genomes = False
    oldies_path = args.genomes_folder
    n_threads = args.threads
    mode = args.mode
    format_genomes=args.format_genomes

    # Genomes ftp selection and download
    ftp_table = gfetcher.ftp_selector(mode=mode, species=guests)

    if host:
        host = ' '.join(host.split(sep='_'))
        host_ftp_table = gfetcher.ftp_selector(mode='single', species=[host])
        ftp_table = ftp_table.append(host_ftp_table, ignore_index=True)

    oldies = gfetcher.fetcher(ftp_table, oldies_path=oldies_path, format_genomes=format_genomes)

    # Database building and indexing
    database = gdatabase.builder(oldies, keep_genomes=keep_genomes, oldies_path=oldies_path)

    idx = galigner.indexer(database, n_threads=n_threads)

    # Alignment and normalization
    alignment = galigner.multi_threaded_aligner(input_folder, idx, mode='basic', n_threads=n_threads)

    norm_alignment = galigner.normalizer(alignment)

    alignment_df = galigner.alignment_to_data_frame(norm_alignment, output_folder=output_folder)

    # Plotting
    barplot.plotter(norm_alignment, alignment_df, host=host, output_folder=output_folder)


if __name__ == '__main__':

    main()

