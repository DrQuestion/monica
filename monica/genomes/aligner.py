import os
import itertools
import pickle
from multiprocessing.dummy import Pool as ThreadPool
from collections import Counter

import pandas as pd
import mappy
from Bio import SeqIO

from .fetcher import GENOMES_PATH
from .database import DATABASE_PATH

BEST_N = 1
INDEX_FILE = os.path.join(DATABASE_PATH, 'index.mmi')

ALIGNMENT_PICKLE_FILENAME = 'alignment.pkl'
ALIGNMENT_PICKLE = os.path.join(os.path.dirname(__file__), ALIGNMENT_PICKLE_FILENAME)

MAPPED_FILES_FOLDER = 'mapped'
UNMAPPED_FILES_FOLDER = 'unmapped'


def indexer(database, n_threads=None, index_file=INDEX_FILE):
    with open(os.path.join(GENOMES_PATH, 'entered_indexer'), 'wb'):
        pass
    index = mappy.Aligner(fn_idx_in=database, preset='map-ont', best_n=BEST_N, n_threads=n_threads, batch_size='250M', fn_idx_out=bytes(index_file, encoding='utf-8'))
    # pickle.dump(index, open(INDEX_PICKLE, 'w'))
    with open(os.path.join(GENOMES_PATH, 'finished_indexing'), 'wb'):
        pass
    return index


def multi_threaded_aligner(query_folder, index, mode=None, overnight=False, n_threads=None,
                           mapped_files_folder=MAPPED_FILES_FOLDER, unmapped_files_folder=UNMAPPED_FILES_FOLDER):

    os.chdir(query_folder)

    samples = [file for file in os.listdir('.') if file.endswith('fastq')]
    samples_name = list(map(lambda sample_name: sample_name.split('.')[0], samples))

    mapped_folder = os.path.join(query_folder, mapped_files_folder)
    unmapped_folder = os.path.join(query_folder, unmapped_files_folder)
    if not os.path.exists(mapped_folder):
        os.mkdir(mapped_folder)
        os.mkdir(unmapped_folder)

    pool = ThreadPool(n_threads)

    results = pool.starmap(aligner, zip(samples, samples_name, itertools.repeat(index), itertools.repeat(mode), itertools.repeat(overnight), itertools.repeat(mapped_folder), itertools.repeat(unmapped_folder)))

    alignment = alignment_update(results)

    return alignment


def aligner(sample, sample_name, index, mode=None, overnight=False, mapped_folder=None, unmapped_folder=None):
    # mode parameter is for testing only
    print(f'{sample}, mode is {mode}\t')
    sample_alignment=dict()
    with open(os.path.join(mapped_folder, sample), 'a') as mapped, \
            open(os.path.join(unmapped_folder, sample), 'a') as unmapped:

        for seq_record in SeqIO.parse(sample, 'fastq'):
            any_hit = 0
            for hit in index.map(str(seq_record.seq)):
                any_hit = 1
                if hit.is_primary:
                    tax_unit = hit.ctg.split(sep=':')[0]
                    if overnight:
                        # tax_unit becomes the genus
                        tax_unit = tax_unit.split(sep='_')[0]
                    accession = hit.ctg.split(sep=':')[1]

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
                            sample_alignment[tax_unit].update({accession: hit.mlen})
                        else:
                            sample_alignment[tax_unit] = Counter({accession: hit.mlen})

            if any_hit:
                SeqIO.write(seq_record, mapped, 'fastq')
            else:
                SeqIO.write(seq_record, unmapped, 'fastq')

    print(f'{sample} done')
    os.remove(sample)
    return sample_alignment, sample_name


def alignment_update(results):
    if ALIGNMENT_PICKLE_FILENAME in os.listdir(os.path.dirname(__file__)):
        alignment=pickle.load(open(ALIGNMENT_PICKLE, 'rb'))
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

    pickle.dump(alignment, open(ALIGNMENT_PICKLE, 'wb'))

    return alignment


def normalizer(alignment):
    genomes_length = pickle.load(open(os.path.join(GENOMES_PATH, 'current_genomes_length.pkl'), 'rb'))
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
