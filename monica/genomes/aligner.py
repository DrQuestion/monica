import os
import gzip
import time
import pickle
from multiprocessing.dummy import Pool as ThreadPool
from collections import Counter

import pandas as pd
import mappy
from Bio import SeqIO

from monica.genomes.fetcher import OLDIES_PATH


BEST_N = 1
N_THREADS=7
IDXFILE = os.path.join(os.path.dirname(__file__), 'index.mmi')

alignment_pickle_filename='alignment.pkl'
alignment_pickle=os.path.join(os.path.dirname(__file__), alignment_pickle_filename)

MAPPED_FILES_FOLDER = 'mapped'
UNMAPPED_FILES_FOLDER = 'unmapped'


def indexer(database, n_threads=N_THREADS, idx_file=False):
    if idx_file:
        index = mappy.Aligner(fn_idx_in=database, preset='map-ont', best_n=BEST_N, n_threads=n_threads, fn_idx_out=IDXFILE)
        return index
    index = mappy.Aligner(fn_idx_in=database, preset='map-ont', best_n=BEST_N, n_threads=n_threads)
    return index


def aligner(query_folder, index, alignment=dict(), mode=None, overnight=False,
            mapped_files_folder=MAPPED_FILES_FOLDER, unmapped_files_folder=UNMAPPED_FILES_FOLDER):
    # mode parameter is for testing only

    # TODO very slow alignment, test with minimap2 cmd line to show how much slow actually
    # TODO look at how multithreading works

    os.chdir(query_folder)

    samples=os.listdir('.')
    samples_name = list(map(lambda sample_name: sample_name.split('.')[0], samples))

    mapped_folder=os.path.join(query_folder, mapped_files_folder)
    unmapped_folder=os.path.join(query_folder, unmapped_files_folder)
    os.mkdir(mapped_folder)
    os.mkdir(unmapped_folder)

    print(samples_name)
    print(alignment)
    for sample, sample_name in zip(samples, samples_name):
        sample_hit_time = 0
        sample_fork_time = 0
        sample_quantification_time = 0
        if sample_name not in alignment:
            alignment[sample_name] = dict()
        with open(os.path.join(mapped_folder, sample), 'w') as mapped, \
                open(os.path.join(unmapped_folder, sample), 'w') as unmapped:

            for seq_record in SeqIO.parse(sample, 'fastq'):
                any_hit = 0
                t_hit=time.time()
                for hit in index.map(str(seq_record.seq)):
                    sample_hit_time += (time.time()-t_hit)
                    any_hit+=1
                    tax_unit = hit.ctg.split(sep=':')[0]
                    if overnight:
                        # taxunit becomes the genus
                        tax_unit=tax_unit.split(sep='_')[0]
                    accession = hit.ctg.split(sep=':')[1]

                    t_quant=time.time()
                    if mode == 'basic':
                        if tax_unit in alignment[sample_name]:
                            alignment[sample_name][tax_unit].update({accession: 1})
                        else:
                            alignment[sample_name][tax_unit] = Counter({accession: 1})
                        sample_quantification_time+=time.time()-t_quant

                    elif mode == 'query_length':
                        if tax_unit in alignment[sample_name]:
                            alignment[sample_name][tax_unit].update({accession: len(seq_record.seq)})
                        else:
                            alignment[sample_name][tax_unit] = Counter({accession: len(seq_record.seq)})

                    elif mode == 'matching':
                        if tax_unit in alignment[sample_name]:
                            alignment[sample_name][tax_unit].update({accession: hit.mlen})
                        else:
                            alignment[sample_name][tax_unit] = Counter({accession: hit.mlen})
                t_bif=time.time()
                if any_hit:
                    SeqIO.write(seq_record, mapped, 'fasta')
                else:
                    SeqIO.write(seq_record, unmapped, 'fasta')
                sample_fork_time+=(time.time()-t_bif)
        print(f'mapped {sample_name}, hitting required {sample_hit_time}, fork required {sample_fork_time}, quantification required {sample_quantification_time}')

    return alignment


def normalizer(alignment):
    genomes_length = pickle.load(open(os.path.join(OLDIES_PATH, 'genomes_length.pkl'), 'rb'))
    for sample in alignment.keys():
        sample_total = 0
        for taxunit, counter in alignment[sample].items():
            for accession, count in counter.items():
                BPB = count/genomes_length[accession]
                sample_total += BPB
                alignment[sample][taxunit][accession] = BPB
        for taxunit, counter in alignment[sample].items():
            for accession, BPB in counter.items():
                BPM = BPB/sample_total
                alignment[sample][taxunit][accession] = BPM
    return alignment


def alignment_to_data_frame(alignment):
    data_frame = pd.concat({k: pd.DataFrame(v).unstack() for k, v in alignment.items()}, axis=1).dropna(how='all')
    return data_frame
