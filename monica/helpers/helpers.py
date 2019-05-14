from re import search, split
from argparse import ArgumentTypeError


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
    units = {'B': 1, 'K': power**2, 'M': power**3, 'G': power**4, 'T': power**5}
    unit = indexing_memory[-1]
    quantity = float(indexing_memory[:-1])
    to_byte = quantity*units[unit]
    return to_byte


if __name__ == '__main__':
    print(human_readable('2.45K'))
