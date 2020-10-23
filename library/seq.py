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
    A={'A'},
    C={'C'},
    G={'G'},
    T={'T'},
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


seq_write_list = ["-"] + [char for char,
                          _ in sorted(seq_read_dict.items(), key=lambda x: x[1])]
