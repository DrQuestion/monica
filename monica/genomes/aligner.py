import os
import gzip
from collections import Counter

import pandas as pd
import mappy
from Bio import SeqIO


BEST_N = 1
IDXFILE = os.path.join(os.path.dirname(__file__), 'index.mmi')

MAPPED_FILES_FOLDER = 'mapped'
UNMAPPED_FILES_FOLDER = 'unmapped'


def indexer(database, idx_file=False):
    if idx_file:
        index = mappy.Aligner(fn_idx_in=database, preset='map-ont', best_n=BEST_N, fn_idx_out=IDXFILE)
        return index
    index = mappy.Aligner(fn_idx_in=database, preset='map-ont', best_n=BEST_N)
    return index


def aligner(query_folder, index, alignment=dict(), genomes_length=dict(), mode=None, overnight=False,
            mapped_files_folder=MAPPED_FILES_FOLDER, unmapped_files_folder=UNMAPPED_FILES_FOLDER):
    # mode parameter is for testing only

    # TODO very slow alignment, test with minimap2 cmd line to show how much slow actually

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
        if sample_name not in alignment:
            alignment[sample_name] = dict()
        with gzip.open(sample, 'rt') as to_map, gzip.open(os.path.join(mapped_folder, sample), 'wt') as mapped, gzip.open(os.path.join(unmapped_folder, sample), 'wt') as unmapped:
            for seq_record in SeqIO.parse(to_map, 'fastq'):
                any_hit = 0
                for hit in index.map(str(seq_record.seq)):
                    any_hit=1
                    tax_unit = hit.ctg.split(sep=':')[0]
                    if overnight:
                        # taxunit becomes the genus
                        tax_unit=tax_unit.split(sep='_')[0]
                    accession = hit.ctg.split(sep=':')[1]
                    if accession not in genomes_length:
                        genomes_length[accession] = hit.ctg_len

                    if mode == 'basic':
                        if tax_unit in alignment[sample_name]:
                            alignment[sample_name][tax_unit].update({accession: 1})
                        else:
                            alignment[sample_name][tax_unit] = Counter({accession: 1})

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

                if any_hit:
                    SeqIO.write(seq_record, mapped, 'fasta')
                else:
                    SeqIO.write(seq_record, unmapped, 'fasta')
        print(f'mapped {sample_name}')

    return alignment, genomes_length


def normalize(alignment, genomes_length):
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
