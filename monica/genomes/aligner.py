import os
import itertools
import pickle
from multiprocessing.dummy import Pool as ThreadPool
from collections import Counter

import pandas as pd
import mappy
from Bio import SeqIO

from .fetcher import GENOMES_PATH
from .database import DATABASES_PATH

BEST_N = 1
INDEXES_PATH = os.path.join(os.path.dirname(__file__), 'indexes')
INDEX_NAME = ['index', '.mmi']

ALIGNMENT_PICKLE_FILENAME = 'alignment.pkl'

MAPPED_FILES_FOLDER = 'mapped'
UNMAPPED_FILES_FOLDER = 'unmapped'
AMBIGUOUS_FILES_FOLDER = 'ambiguous'
FOCUS_FILES_FOLDER = 'focus'


def indexer(databases, indexes_path=INDEXES_PATH):
    if not os.path.exists(indexes_path):
        os.mkdir(indexes_path)
    else:
        for index in os.listdir(indexes_path):
            os.remove(os.path.join(indexes_path, index))
    indexes = []

    with open(os.path.join(GENOMES_PATH, 'entered_indexer'), 'wb'):
        pass
    for database in os.listdir(databases):
        if database.endswith('.fna.gz'):
            number = os.path.split(database)[1][8:-7]
            index = mappy.Aligner(fn_idx_in=os.path.join(DATABASES_PATH, database), preset='map-ont', best_n=BEST_N,
                                  fn_idx_out=bytes(os.path.join(INDEXES_PATH, str(number).join(INDEX_NAME)),
                                                   encoding='utf-8'))
            if not index:
                raise Exception('Index building failed')
            indexes.append(index)
    with open(os.path.join(GENOMES_PATH, 'finished_indexing'), 'wb'):
        pass
    return indexes


def indexes_opener(indexes_path=INDEXES_PATH):
    indexes = []
    for file in os.listdir(indexes_path):
        if file.endswith('.mmi'):
            index = mappy.Aligner(fn_idx_in=os.path.join(indexes_path, file))
            indexes.append(index)
            if not index:
                raise Exception('Damaged or empty index')
    return indexes


def multi_threaded_aligner(query_folder, indexes, mode=None, overnight=False, n_threads=None, focus_species=[],
                           output_folder=None, mapped_files_folder=MAPPED_FILES_FOLDER,
                           unmapped_files_folder=UNMAPPED_FILES_FOLDER, ambiguous_files_folder=AMBIGUOUS_FILES_FOLDER,
                           focus_file_folder=FOCUS_FILES_FOLDER):

    os.chdir(query_folder)

    samples = [file for file in os.listdir('.') if file.endswith('fastq')]
    samples_name = list(map(lambda sample_name: sample_name.split('.')[0], samples))

    mapped_folder = os.path.join(query_folder, mapped_files_folder)
    unmapped_folder = os.path.join(query_folder, unmapped_files_folder)
    ambiguous_folder = os.path.join(query_folder, ambiguous_files_folder)
    focus_folder = os.path.join(query_folder, focus_file_folder)
    if not os.path.exists(mapped_folder):
        os.mkdir(mapped_folder)
        os.mkdir(unmapped_folder)
        os.mkdir(ambiguous_folder)
        if focus_species:
            os.mkdir(focus_folder)

    pool = ThreadPool(n_threads)

    results = pool.starmap(aligner, zip(samples, samples_name, itertools.repeat(indexes), itertools.repeat(mode),
                                        itertools.repeat(overnight), itertools.repeat(focus_species),
                                        itertools.repeat(mapped_folder), itertools.repeat(unmapped_folder),
                                        itertools.repeat(ambiguous_folder), itertools.repeat(focus_folder)))

    alignment = alignment_update(results, output_folder)

    return alignment


def aligner(sample, sample_name, indexes, mode=None, overnight=False, focus_species=[], mapped_folder=None,
            unmapped_folder=None, ambiguous_folder=None, focus_folder=None):
    # mode parameter is for testing only
    print(f'{sample}, mode is {mode}\t')
    sample_alignment = dict()

    if focus_species:
        focus = open(os.path.join(focus_folder, sample), 'a')

    with open(os.path.join(mapped_folder, sample), 'a') as mapped, \
            open(os.path.join(unmapped_folder, sample), 'a') as unmapped, \
            open(os.path.join(ambiguous_folder, sample), 'a') as ambiguous:

        for seq_record in SeqIO.parse(sample, 'fastq'):
            hits = []
            for index in indexes:
                for hit in index.map(str(seq_record.seq)):
                    if hit.is_primary:
                        hits.append(hit)
            if hits:
                if len(hits) == 1:
                    best = hits[0]
                else:
                    best = best_hit(hits)
                    if not best:
                        SeqIO.write(seq_record, ambiguous, 'fastq')
                        continue
                SeqIO.write(seq_record, mapped, 'fastq')
                tax_unit = best.ctg.split(sep=':')[0]
                if tax_unit in focus_species:
                    SeqIO.write(seq_record, focus, 'fastq')
                if overnight:
                    # tax_unit becomes the genus
                    tax_unit = tax_unit.split(sep='_')[0]
                accession = best.ctg.split(sep=':')[1]

                # TODO: Implement a way to keep them all and plot them alternatively at the end

                if mode == 'basic':
                    if tax_unit in sample_alignment:
                        sample_alignment[tax_unit].update({accession: 1})
                    else:
                        sample_alignment[tax_unit] = Counter({accession: 1})

                elif mode == 'query_length':
                    if tax_unit in sample_alignment:
                        sample_alignment[tax_unit].update({accession: len(seq_record.seq)})
                    else:
                        sample_alignment[tax_unit] = Counter({accession: len(seq_record.seq)})

                elif mode == 'matching':
                    if tax_unit in sample_alignment:
                        sample_alignment[tax_unit].update({accession: best.mlen})
                    else:
                        sample_alignment[tax_unit] = Counter({accession: best.mlen})
            else:
                SeqIO.write(seq_record, unmapped, 'fastq')

    if focus_species:
        focus.close()
    print(f'{sample} done')
    os.remove(sample)
    return sample_alignment, sample_name


def alignment_update(results, output_folder):
    alignment_pickle = os.path.join(output_folder, ALIGNMENT_PICKLE_FILENAME)
    if ALIGNMENT_PICKLE_FILENAME in os.listdir(output_folder):
        alignment = pickle.load(open(alignment_pickle, 'rb'))
        for alignment_sample, sample_name in results:
            if sample_name in alignment:
                for tax_unit, counter in alignment_sample.items():
                        if tax_unit in alignment[sample_name]:
                            alignment[sample_name][tax_unit].update(counter)
                        else:
                            alignment[sample_name][tax_unit] = counter
            else:
                alignment[sample_name] = alignment_sample
    else:
        alignment = dict()
        for alignment_sample, sample_name in results:
            alignment[sample_name] = alignment_sample

    pickle.dump(alignment, open(alignment_pickle, 'wb'))

    return alignment


def normalizer(alignment, genomes_length=None, databases_path=None):
    if not genomes_length:
        genomes_length = pickle.load(open(os.path.join(databases_path, 'current_genomes_length.pkl'), 'rb'))
    for sample in alignment.keys():
        sample_total = 0
        for tax_unit, counter in alignment[sample].items():
            for accession, count in counter.items():
                BPB = count/genomes_length[accession]
                sample_total += BPB
                alignment[sample][tax_unit][accession] = BPB
        for tax_unit, counter in alignment[sample].items():
            for accession, BPB in counter.items():
                BPM = BPB/sample_total
                alignment[sample][tax_unit][accession] = BPM
    return alignment


def alignment_to_data_frame(alignment, output_folder=None, filename='monica.dataframe'):
    data_frame = pd.concat({k: pd.DataFrame(v).unstack() for k, v in alignment.items()}, axis=1).dropna(how='all')
    pd.DataFrame.to_csv(data_frame, os.path.join(output_folder, filename))
    return data_frame


def best_hit(hits):
    best = (None, float("inf"))
    distance = 0
    for hit in hits:
        inverse_identity = float(hit.NM)/hit.mlen
        if inverse_identity <= best[1]:
            distance = best[1] - inverse_identity
            best = (hit, inverse_identity)
    if not distance:
        return 0
    else:
        return best[0]
