import os
from collections import Counter

import pandas as pd
import mappy


IDXFILE=os.path.join(os.path.dirname(__file__), 'index.mmi')


def indexer (database, idx_file=False):
    if idx_file:
        index=mappy.Aligner(fn_idx_in=database, preset='map-ont', fn_idx_out=IDXFILE)
        return index
    index = mappy.Aligner(fn_idx_in=database, preset='map-ont')
    return index


def aligner (query_folder, index, alignment=dict(), genomes_length=dict(), mode=None, overnight=False):
    # mode is for testing only

    os.chdir(query_folder)
    samples=os.listdir('.')

    for sample in samples:
        alignment[sample]=dict()
        for name, seq, qual in mappy.fastx_read(sample):
            for hit in index.map(seq):
                # When do I cut secondary alignments?
                tax_unit = hit.ctg.split(sep=':')[0]
                if overnight:
                    # taxunit becomes the genus
                    tax_unit=tax_unit.split(sep='_')[0]
                accession = hit.ctg.split(sep=':')[1]
                if not accession in genomes_length:
                    genomes_length[accession] = hit.ctg_len

                if mode=='basic':
                    if tax_unit in alignment[sample]:
                        alignment[sample][tax_unit].update({accession: 1})
                    else:
                        alignment[sample][tax_unit] = Counter({accession: 1})

                elif mode=='query_length':
                    if tax_unit in alignment[sample]:
                        alignment[sample][tax_unit].update({accession: len(seq)})
                    else:
                        alignment[sample][tax_unit] = Counter({accession: len(seq)})

                elif mode=='matching':
                    if tax_unit in alignment[sample]:
                        alignment[sample][tax_unit].update({accession: hit.mlen})
                    else:
                        alignment[sample][tax_unit] = Counter({accession: hit.mlen})

    return alignment, genomes_length


def normalize(alignment, genomes_length):
    for sample in alignment.keys():
        sample_total=0
        for taxunit, counter in alignment[sample].items():
            for accession, count in counter.items():
                BPB=count/genomes_length[accession]
                sample_total+=BPB
                alignment[sample][taxunit][accession]=BPB
        for taxunit, counter in alignment[sample].items():
            for accession, BPB in counter.items():
                BPM=BPB/sample_total
                alignment[sample][taxunit][accession]=BPM
    return alignment


def allignment_to_data_frame (alignment):
    data_frame = pd.concat({k: pd.DataFrame(v).unstack() for k, v in alignment.items()}, axis=1)
    return data_frame


if __name__=='__main__':
    db='/home/drq/Desktop/temp/genome.fna'
    query='/home/drq/Desktop/temp/query'
    al=aligner(query, db)
    print(al)
