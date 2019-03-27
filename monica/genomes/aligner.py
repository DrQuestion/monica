from collections import Counter
import mappy as mp
db='/home/drq/Desktop/genome.fna'
query='/home/drq/Desktop/query.fa'

if __name__=='__main__':
    g_lengths=dict()
    tax_units=dict()
    a=mp.Aligner(db)
    with open (query, 'r') as f:
        print(f.read())
    for name, seq, qual in mp.fastx_read(query):
        print(len(seq))
        for hit in a.map(seq):
            tax_unit=hit.ctg.split(sep=':')[0]
            accession=hit.ctg.split(sep=':')[1]
            if tax_unit in tax_units:
                tax_units[tax_unit].update({accession:1})
            else:
                tax_units[tax_unit]=Counter({accession:1})
                g_lengths[accession]=hit.ctg_len
            print(tax_units)
            print(g_lengths)