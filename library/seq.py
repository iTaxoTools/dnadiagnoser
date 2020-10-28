from typing import Union, Iterator, Tuple
import numpy as np
from functools import reduce
seq_read_dict = dict(
    A=1,
    C=2,
    G=4,
    T=8,
    R=5,
    Y=10,
    S=6,
    W=9,
    K=12,
    M=3,
    B=14,
    D=13,
    H=11,
    V=7,
    N=15
)

seq_set_dict = dict(
    A=frozenset('A'),
    C=frozenset('C'),
    G=frozenset('G'),
    T=frozenset('T'),
    R=frozenset("AG"),
    Y=frozenset("CT"),
    S=frozenset("GC"),
    W=frozenset("AT"),
    K=frozenset("GT"),
    M=frozenset("AC"),
    B=frozenset("CGT"),
    D=frozenset("AGT"),
    H=frozenset("ACT"),
    V=frozenset("ACG"),
    N=frozenset("ACGT"),
)

assert(sorted(seq_read_dict.keys()) == sorted(seq_set_dict.keys()))

assert(all(seq_read_dict[key] == reduce(lambda acc, x: acc |
                                        seq_read_dict[x], seq_set_dict[key], 0) for key in seq_read_dict.keys()))


seq_write_tuple = ("-",) + tuple(char for char,
                                 _ in sorted(seq_read_dict.items(), key=lambda x: x[1]))


class Seq:
    """
    Store a sequence in a compact way as np.array of bytes

    Each byte encodes a set of nucleotides into the 4 lower bits with correspondences:
    A <~> 0b0001
    C <~> 0b0010
    G <~> 0b0100
    T <~> 0b1000
    """

    def __init__(self, data: np.array) -> None:
        self.data = data

    @classmethod
    def from_str(cls, sequence: str) -> 'Seq':
        seq = cls(np.empty(len(sequence), dtype="byte"))
        for i, el in enumerate(map(seq_read_dict.get, sequence)):
            seq.data[i] = el
        return seq

    def __or__(self, other: 'Seq') -> 'Seq':
        return Seq(self.data | other.data)

    def __str__(self) -> str:
        return "".join(map(seq_write_tuple.__getitem__, self.data))

    def __iter__(self) -> Iterator:
        return self.data.__iter__()


def differences(seq1: Seq, seq2: Seq) -> Iterator[Tuple[int, int, int]]:
    """
    Returns an iterator over incompatible nucleotides of arguments with indices

    Iterator element: (index, nucleotide1, nucleotide2)
    """
    return ((i, nuc1, nuc2) for i, (nuc1, nuc2) in enumerate(zip(seq1, seq2)) if not nuc1 & nuc2 and nuc1 and nuc2)
