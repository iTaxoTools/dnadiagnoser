#!/usr/bin/env python3

from library.seq import Seq, seq_write_tuple, differences, Dict, List, Tuple
import sys
import pandas as pd
import numpy as np
from functools import reduce
import itertools
import os.path
from library.gui_utils import *
import tkinter.messagebox


with open(os.path.join('data', 'dnadiagnoser_homo_sapiens_reference.fas')) as file:
    references: Dict[str, Seq] = {}
    for line in file:
        name, _, sequence = line.partition('\t')
        if name and sequence:
            references[name] = Seq.from_str(sequence)


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


def process_files(infile: str, outfile: str, reference_name: str) -> None:
    if not infile:
        raise ValueError('Input file is not given')
    if not outfile:
        raise ValueError('Output file is not given')
    table = pd.read_csv(infile, delimiter='\t', index_col='specimen_voucher', converters={
                        'sequence': Seq.from_str})
    reference_sequence = references[reference_name]
    alignment_displays = table['sequence'].apply(
        lambda seq: seq.align(reference_sequence))
    table = table.groupby('species')['sequence'].agg(combine_sequences)
    with open(outfile, mode='w') as output:
        print("Alignments:", file=output)
        for specimen, alignment_str in alignment_displays.items():
            print(specimen, file=output)
            print(alignment_str, file=output)
        output.write("\n")
        print("species 1\tspecies 2\treplacements\tinsertions 1\tinsertions 2", file=output)
        for species1, species2 in itertools.product(table.index, table.index):
            if species1 == species2:
                continue
            repl, ins1, ins2 = differences(table[species1], table[species2])
            print(species1, species2, show_differences(
                repl, ins1, ins2), sep='\t', file=output)


def launch_gui() -> None:
    root = tk.Tk()
    root.rowconfigure(0, weight=1)
    root.columnconfigure(0, weight=1)
    root.columnconfigure(2, weight=1)

    inputchooser = FileChooser(root, label="Input file", mode="open")
    outputchooser = FileChooser(root, label="Output file", mode="save")

    reference_cmb = LabeledCombobox(
        root, label="Reference sequence", values=list(references.keys()), readonly=True)

    def process() -> None:
        try:
            process_files(inputchooser.file_var.get(),
                          outputchooser.file_var.get(), reference_cmb.var.get())
        except Exception as ex:
            tkinter.messagebox.showerror("Error", str(ex))

    process_btn = ttk.Button(root, text="Process", command=process)

    inputchooser.grid(row=0, column=0, sticky="nsew")
    outputchooser.grid(row=0, column=2, sticky="nsew")

    reference_cmb.grid(row=1, column=1)
    process_btn.grid(row=2, column=1)

    root.mainloop()


def main() -> None:
    if len(sys.argv) >= 3:
        input = sys.argv[1]
        output = sys.argv[2]
        try:
            reference_name = sys.argv[3]
        except IndexError:
            reference_name = "Homo_sapiens_COI"
        process_files(input, output, reference_name)
    else:
        launch_gui()


if __name__ == "__main__":
    main()
