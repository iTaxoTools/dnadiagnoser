#!/usr/bin/env python3

import io
import itertools
import os.path
import sys
import tkinter.messagebox
import warnings
from functools import reduce
from typing import Dict, List, Optional, TextIO, Tuple

import numpy as np
import pandas as pd

from library.gui_utils import *
from library.seq import Seq, differences, seq_write_tuple

with open(os.path.join('data', 'dnadiagnoser_homo_sapiens_reference.fas')) as file:
    references: Dict[str, Seq] = {}
    for line in file:
        name, _, sequence = line.partition('\t')
        if name and sequence:
            references[name] = Seq.from_str(sequence)

typos = dict(
    specimen_voucher='specimenid',
    specimen_id='specimenid',
    sequences='sequence'
)


def combine_sequences(series: pd.Series) -> Seq:
    return reduce(lambda s1, s2: s1 | s2, series)


def and_join(words: List[str]) -> str:
    if not words:
        return ""
    elif len(words) == 1:
        return words[0]
    else:
        return ", ".join(words[:-1]) + " and " + words[-1]


def show_differences(repl: List[Tuple[int, int, int]], ins1: Dict[int, np.array], ins2: Dict[int, np.array]) -> str:
    replacements = ", ".join(
        f"{i} ({seq_write_tuple[nuc1]} vs. {seq_write_tuple[nuc2]})" for i, nuc1, nuc2 in repl)
    insertions1 = ", ".join(
        f"at {i} {''.join(seq_write_tuple[nuc] for nuc in frag)}" for i, frag in ins1.items())
    insertions2 = ", ".join(
        f"at {i} {''.join(seq_write_tuple[nuc] for nuc in frag)}" for i, frag in ins2.items())
    return "\t".join((replacements, insertions1, insertions2))


def textual_differences(repl: List[Tuple[int, int, int]], ins1: Dict[int, np.array], ins2: Dict[int, np.array]) -> str:
    repls_as_str = [(i, f"({seq_write_tuple[nuc1]} vs. {seq_write_tuple[nuc2]})")
                    for i, nuc1, nuc2 in repl]
    ins_as_str = [
        (i, f"(insertion {''.join(seq_write_tuple[nuc] for nuc in frag)})") for i, frag in ins1.items()]
    del_as_str = [
        (i, f"(deletion {''.join(seq_write_tuple[nuc] for nuc in frag)})") for i, frag in ins1.items()]
    listed_difference = [f"{i} {desc}" for i, desc in sorted(
        repls_as_str + ins_as_str + del_as_str, key=lambda x: x[0])]
    return and_join(listed_difference)


def show_diag_differences(repl: List[Tuple[int, int, int]], ins1: Dict[int, np.array], ins2: Dict[int, np.array]) -> str:
    repls_as_str = [(i, f"({seq_write_tuple[nuc1]})")
                    for i, nuc1, nuc2 in repl]
    ins_as_str = [
        (i, f"(insertion {''.join(seq_write_tuple[nuc] for nuc in frag)})") for i, frag in ins1.items()]
    del_as_str = [
        (i, f"(deletion {''.join(seq_write_tuple[nuc] for nuc in frag)})") for i, frag in ins1.items()]
    listed_difference = [f"{i} {desc}" for i, desc in sorted(
        repls_as_str + ins_as_str + del_as_str, key=lambda x: x[0])]
    return ", ".join(listed_difference)


def diag_textual_differences(repl: List[Tuple[int, int, int]], ins1: Dict[int, np.array], ins2: Dict[int, np.array]) -> str:
    repls_as_str = [(i, f"having a {seq_write_tuple[nuc1]}")
                    for i, nuc1, nuc2 in repl]
    ins_as_str = [
        (i, f"having an insertion {''.join(seq_write_tuple[nuc] for nuc in frag)}") for i, frag in ins1.items()]
    del_as_str = [
        (i, f"having a deletion {''.join(seq_write_tuple[nuc] for nuc in frag)}") for i, frag in ins1.items()]
    listed_difference = [f"{desc} at position {i}" for i, desc in sorted(
        repls_as_str + ins_as_str + del_as_str, key=lambda x: x[0])]
    return and_join(listed_difference)


class DnaProcessor():

    def __init__(self) -> None:
        self.infile: Optional[str] = None
        self.table: Optional[pd.DataFrame] = None

    def load_table(self, infile: str) -> None:
        table = pd.read_csv(infile, delimiter='\t').rename(
            columns=str.casefold).rename(columns=typos)
        if 'sequence' not in table.columns:
            raise ValueError("'sequences' or 'sequence' column is missing")
        if 'specimenid' not in table.columns:
            warnings.warn("Specimen IDs are not detected")
        else:
            table.set_index('specimenid', inplace=True)
        if len(table.columns) < 2:
            raise ValueError("'species' or another column need to be present")
        table['sequence'] = table['sequence'].apply(Seq.from_str)
        self.table = table
        self.infile = infile

    def process_files(self, infile: str, outfile: str, reference_name: str, column: str, selection: List[str]) -> None:
        if not infile:
            raise ValueError('Input file is not given')
        if not outfile:
            raise ValueError('Output file is not given')
        if self.infile != infile or self.table is None:
            self.load_table(infile)
            column = "species"
            selection = []
        assert(self.table is not None)
        reference_sequence = references[reference_name]
        alignment_displays = self.table['sequence'].apply(
            lambda seq: seq.align(reference_sequence))
        table = self.table.groupby(
            column)['sequence'].agg(combine_sequences)
        if selection:
            table = table.filter(items=selection).sort_index()
        with open(outfile, mode='w') as output:
            print("Alignments:", file=output)
            for specimen, alignment_str in alignment_displays.items():
                if not selection or self.table[column][specimen] in selection:
                    print(specimen, file=output)
                    print(alignment_str, file=output)
            output.write("\n")
            self.report(output, table, column, reference_name)

    def report(self, output: TextIO, table: pd.Series, column: str, reference_name: str) -> None:
        matrixOutput = io.StringIO("")
        tableOutput = io.StringIO("")
        textOutput = io.StringIO("")
        print("", *table.index, sep='\t', file=matrixOutput)
        print(
            f"{column} 1\t{column} 2\treplacements\tinsertions 1\tinsertions 2", file=tableOutput)
        num_species = len(table.index)
        for species1 in table.index:
            matrixOutput.write(species1)
            textOutput.write(
                f"Using nucleotide positions in the {reference_name} sequence as a reference, {species1} differs ")
            textFragments = []
            for species2 in table.index:
                repl, ins1, ins2 = differences(
                    table[species1], table[species2])
                insertions_num = sum(
                    map(len, ins1.values())) + sum(map(len, ins2.values()))
                difference_num = len(repl) + insertions_num
                if species1 >= species2:
                    matrixOutput.write(f"\t{difference_num}")
                if species1 != species2:
                    print(species1, species2, show_differences(
                        repl, ins1, ins2), sep='\t', file=tableOutput)
                if difference_num > 0:
                    textFragments.append(
                        f"from {species2} in nucleotide {'position' if difference_num == 1 else 'positions'} " +
                        textual_differences(repl, ins1, ins2))
            textOutput.write(
                and_join(textFragments))
            textOutput.write("\n\n")
            matrixOutput.write("\n")
        print(matrixOutput.getvalue(), file=output)
        print(tableOutput.getvalue(), file=output)
        print(textOutput.getvalue(), file=output)

        diag_tableOutput = io.StringIO("")
        diag_textOutput = io.StringIO("")
        print(f"{column}\tUnique diagnostic differences", file=diag_tableOutput)
        for species1 in table.index:
            other_species_seq: Seq = combine_sequences(table.drop(
                labels=species1))
            repl, ins1, ins2 = differences(table[species1], other_species_seq)
            if (text := show_diag_differences(repl, ins1, ins2)):
                print(species1, text, sep='\t', file=diag_tableOutput)
            else:
                print(species1, "None", sep='\t', file=diag_tableOutput)
            if (text := diag_textual_differences(repl, ins1, ins2)):
                diag_textOutput.write(
                    f"{species1} differs from all other {column} in the dataset by ")
                diag_textOutput.write(text)
                diag_textOutput.write(
                    f" of the {reference_name} reference sequence.\n")
            else:
                print(
                    f"{species1} has no unique diagnostic differences in comparison to the other {column} in the data set.", file=diag_textOutput)
        print(diag_tableOutput.getvalue(), file=output)
        print(diag_textOutput.getvalue(), file=output)

    def choices(self) -> Dict[str, List[str]]:
        """
        Returns the columns that contain values for grouping
        """
        assert(self.table is not None)
        return {column_name: list(self.table[column_name].unique()) for column_name in self.table.columns if column_name != "sequence"}


def launch_gui() -> None:
    root = tk.Tk()
    root.rowconfigure(0, weight=1)
    root.columnconfigure(0, weight=1)
    root.columnconfigure(2, weight=1)

    inputchooser = FileChooser(root, label="Input file", mode="open")
    outputchooser = FileChooser(root, label="Output file", mode="save")

    reference_cmb = LabeledCombobox(
        root, label="Reference sequence", values=list(references.keys()), readonly=True)

    processor = DnaProcessor()

    column_selector = ColumnSelector(root)

    def load() -> None:
        try:
            if (infile := inputchooser.file_var.get()):
                processor.load_table(infile)
                column_selector.set_columns(processor.choices())
        except Exception as ex:
            tkinter.messagebox.showerror("Error", str(ex))

    def process() -> None:
        try:
            with warnings.catch_warnings(record=True) as warns:
                if (selection := column_selector.selection()):
                    column, selected = selection
                else:
                    column, selected = ('species', [])
                processor.process_files(inputchooser.file_var.get(),
                                        outputchooser.file_var.get(), reference_cmb.var.get(), column, selected)
                for w in warns:
                    tkinter.messagebox.showwarning("Warning", str(w.message))
        except Exception as ex:
            tkinter.messagebox.showerror("Error", str(ex))
        else:
            tkinter.messagebox.showinfo("Done", "Analysis is complete")

    process_btn = ttk.Button(root, text="Process", command=process)
    load_btn = ttk.Button(root, text="Load", command=load)

    inputchooser.grid(row=0, column=0, sticky="nsew")
    outputchooser.grid(row=0, column=2, sticky="nsew")

    reference_cmb.grid(row=1, column=1)
    process_btn.grid(row=2, column=1)
    load_btn.grid(row=2, column=0)
    column_selector.grid(row=3, column=0, sticky="nsew")

    root.mainloop()


def main() -> None:
    if len(sys.argv) >= 3:
        input = sys.argv[1]
        output = sys.argv[2]
        try:
            reference_name = sys.argv[3]
        except IndexError:
            reference_name = "Homo_sapiens_COI"
        processor = DnaProcessor()
        processor.process_files(input, output, reference_name, "species", [])
    else:
        launch_gui()


if __name__ == "__main__":
    main()
