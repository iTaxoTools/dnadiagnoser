from typing import Union, Iterator, Tuple, Dict, List
import numpy as np
from functools import reduce
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
import itertools
import os

with open(os.path.join('data', 'scores.tab')) as scores_file:
    scores_dict = {}
    for line in scores_file:
        score_name, _, val = line.partition('\t')
        try:
            scores_dict[score_name] = int(val)
        except ValueError as ex:
            raise ValueError(
                f"The value for '{score_name}' in data/scores.tab is not a number") from ex

try:
    GAP_PENALTY = scores_dict['gap penalty']
    GAP_EXTEND_PENALTY = scores_dict['gap extend penalty']
    END_GAP_PENALTY = scores_dict['end gap penalty']
    END_GAP_EXTEND_PENALTY = scores_dict['end gap extend penalty']
    MATCH_SCORE = scores_dict['match score']
    MISMATCH_SCORE = scores_dict['mismatch score']
except KeyError as ex:
    raise ValueError(f"'{ex.args[0]}' is missing in data/scores.tab") from ex

score_matrix = np.empty((16, 16))

for i, j in itertools.product(range(0, 16), range(0, 16)):
    score_matrix[i, j] = MATCH_SCORE if i & j else MISMATCH_SCORE

aligner = PairwiseAligner(substitution_matrix=score_matrix, end_open_gap_score=END_GAP_PENALTY,
                          end_extend_gap_score=END_GAP_EXTEND_PENALTY, internal_open_gap_score=GAP_PENALTY, internal_extend_gap_score=GAP_EXTEND_PENALTY)

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


def np_or_maxlen(arr1: np.array, arr2: np.array) -> np.array:
    if len(arr1) < len(arr2):
        arr1.resize(len(arr2))
    elif len(arr2) < len(arr1):
        arr2.resize(len(arr1))
    return arr1 | arr2


def merge_insertions(ins1: Dict[int, np.array], ins2: Dict[int, np.array]) -> Dict[int, np.array]:
    result = {**ins1, **ins2}
    result.update({key: np_or_maxlen(ins1[key], ins2[key])
                   for key in ins1.keys() & ins2.keys()})
    return result


class Seq:
    """
    Store a sequence in a compact way as np.array of bytes

    Each byte encodes a set of nucleotides into the 4 lower bits with correspondences:
    A <~> 0b0001
    C <~> 0b0010
    G <~> 0b0100
    T <~> 0b1000
    """

    def __init__(self, data: np.array, insertions: Dict[int, np.array]) -> None:
        self.data = data
        self.insertions = insertions

    @classmethod
    def from_str(cls, sequence: str) -> 'Seq':
        sequence = sequence.strip("-Nn?\n\t ").upper()
        seq = cls(np.empty(len(sequence), dtype="int32"), {})
        for i, el in enumerate(map(seq_read_dict.get, sequence)):
            try:
                seq.data[i] = el
            except TypeError as ex:
                raise ValueError(
                    f"Unexpected nucleotide: {sequence[i]}") from ex
        return seq

    def __or__(self, other: 'Seq') -> 'Seq':
        return Seq(self.data | other.data, merge_insertions(self.insertions, other.insertions))

    # def __str__(self) -> str:
    #     return "".join(map(seq_write_tuple.__getitem__, self.data))

    def __iter__(self) -> Iterator:
        return self.data.__iter__()

    def align(self, ref: 'Seq') -> str:
        alignment = aligner.align(ref.data, self.data)[0]
        aligned = alignment.aligned
        alignment.target = "".join(
            map(seq_write_tuple.__getitem__, alignment.target))
        alignment.query = "".join(
            map(seq_write_tuple.__getitem__, alignment.query))
        aligned_data = np.zeros_like(ref.data)
        self.insertions = {}
        prev_self_end = 0
        prev_ref_end = 0
        for (ref_start, ref_end), (self_start, self_end) in zip(*aligned):
            aligned_data[ref_start:ref_end] = self.data[self_start:self_end]
            if self_start > prev_self_end:
                self.insertions[prev_ref_end] = self.data[prev_self_end: self_start]
            prev_self_end = self_end
            preve_ref_end = ref_end
        else:
            if prev_self_end < len(self.data):
                self.insertions[prev_ref_end] = self.data[prev_self_end:]
        self.data = aligned_data
        return format(alignment)


def differences(seq1: Seq, seq2: Seq) -> Tuple[List[Tuple[int, int, int]], Dict[int, np.array], Dict[int, np.array]]:
    """
    Returns a list of (index, nucleotide1, nucleotide2) of different nucleotides.

    Additionally returns a pair of dictionaries for different insertions for each sequence
    """
    replacements = [(i, nuc1, nuc2) for i, (nuc1, nuc2) in enumerate(
        zip(seq1, seq2)) if not nuc1 & nuc2 and (nuc1 or nuc2)]
    ins1 = seq1.insertions
    ins2 = seq2.insertions
    common_insertions = {key: (ins1[key], ins2[key])
                         for key in ins1.keys() & ins2.keys()}
    ins1 = {key: ins1[key] for key in ins1 if key not in ins2}
    ins1.update({key: val[0] for key, val in common_insertions.items()})
    ins2 = {key: ins2[key] for key in ins2 if key not in ins1}
    ins2.update({key: val[0] for key, val in common_insertions.items()})

    return replacements, ins1, ins2
