import os
from re import search
from argparse import ArgumentTypeError

from genomes.aligner import ALIGNMENT_PICKLE


def human_readable(indexing_memory):
    if search('[BKMGT]$', indexing_memory):
        pass
    else:
        msg = 'Wrong memory unit specified, it must end with one among B|K|M|G|T'
        raise ArgumentTypeError(msg)
    if search('^\d', indexing_memory):
        return human_readable_to_byte(indexing_memory)
    else:
        msg = 'Should also put a number in front of memory unit..'
        raise ArgumentTypeError(msg)


def human_readable_to_byte(indexing_memory):
    power = 2**10
    units = {'B': 1, 'K': power**1, 'M': power**2, 'G': power**3, 'T': power**4}
    unit = indexing_memory[-1]
    quantity = float(indexing_memory[:-1])
    to_byte = quantity*units[unit]
    return to_byte


def initializer():
    if os.path.exists(ALIGNMENT_PICKLE):
        os.remove(ALIGNMENT_PICKLE)


if __name__ == '__main__':
    print(human_readable('2.45K'))
