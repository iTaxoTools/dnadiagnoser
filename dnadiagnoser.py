#!/usr/bin/env python3

from library.seq import Seq, seq_write_tuple, differences, Dict, List, Tuple
import sys
import pandas as pd
import numpy as np
from functools import reduce
import itertools
import os.path


with open(os.path.join('data', 'dnadiagnoser_homo_sapiens_reference.fas')) as file:
    name, seq = file.readline().split()
    reference_sequence = Seq.from_str(seq)


def combine_sequences(series: pd.Series) -> Seq:
    return reduce(lambda s1, s2: s1 | s2, series)


def show_differences(repl: List[Tuple[int, int, int]], ins1: Dict[int, np.array], ins2: Dict[int, np.array]) -> str:
    replacements = ", ".join(
        f"{i} ({seq_write_tuple[nuc1]} vs. {seq_write_tuple[nuc2]})" for i, nuc1, nuc2 in repl)
    insertions1 = ", ".join(
        f"at {i} {''.join(seq_write_tuple[nuc] for nuc in frag)}" for i, frag in ins1.items())
    insertions2 = ", ".join(
        f"at {i} {''.join(seq_write_tuple[nuc] for nuc in frag)}" for i, frag in ins2.items())
    return "\t".join((replacements, insertions1, insertions2))


def main() -> None:
    input = sys.argv[1]
    outfile = sys.argv[2]

    table = pd.read_csv(input, delimiter='\t', index_col='specimen_voucher', converters={
                        'sequence': Seq.from_str})
    table['sequence'].apply(
        lambda seq: seq.align(reference_sequence))
    table = table.groupby('species')['sequence'].agg(combine_sequences)
    with open(outfile, mode='w') as output:
        print("species 1\tspecies 2\treplacements\tinsertions 1\tinsertions 2", file=output)
        for species1, species2 in itertools.product(table.index, table.index):
            if species1 == species2:
                continue
            repl, ins1, ins2 = differences(table[species1], table[species2])
            print(species1, species2, show_differences(
                repl, ins1, ins2), sep='\t', file=output)


if __name__ == "__main__":
    main()
