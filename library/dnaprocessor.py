import os
import io
import tkinter as tk
from typing import Dict, List, Optional, TextIO, Callable, Tuple
import warnings
from functools import reduce

import pandas as pd
import numpy as np

from library.seq import Seq, differences, seq_write_tuple

with open(os.path.join('data', 'reference_sequences.tab')) as file:
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
    return reduce(Seq.__or__, series)


def and_join(words: List[str]) -> str:
    if not words:
        return ""
    elif len(words) == 1:
        return words[0]
    else:
        return ", ".join(words[:-1]) + " and " + words[-1]


def show_differences(repl: List[Tuple[int, int, int]], ins1: Dict[int, np.array], ins2: Dict[int, np.array], translation: Callable[[int], str] = str) -> str:
    replacements = ", ".join(
        f"{translation(i)} ({seq_write_tuple[nuc1]} vs. {seq_write_tuple[nuc2]})" for i, nuc1, nuc2 in repl)
    insertions1 = ", ".join(
        f"at {translation(i)} {''.join(seq_write_tuple[nuc] for nuc in frag)}" for i, frag in ins1.items())
    insertions2 = ", ".join(
        f"at {translation(i)} {''.join(seq_write_tuple[nuc] for nuc in frag)}" for i, frag in ins2.items())
    return "\t".join((replacements, insertions1, insertions2))


def textual_differences(repl: List[Tuple[int, int, int]], ins1: Dict[int, np.array], ins2: Dict[int, np.array], translation: Callable[[int], str] = str) -> str:
    repls_as_str = [(i, f"({seq_write_tuple[nuc1]} vs. {seq_write_tuple[nuc2]})")
                    for i, nuc1, nuc2 in repl]
    ins_as_str = [
        (i, f"(insertion {''.join(seq_write_tuple[nuc] for nuc in frag)})") for i, frag in ins1.items()]
    del_as_str = [
        (i, f"(deletion {''.join(seq_write_tuple[nuc] for nuc in frag)})") for i, frag in ins2.items()]
    listed_difference = [f"{translation(i)} {desc}" for i, desc in sorted(
        repls_as_str + ins_as_str + del_as_str, key=lambda x: x[0])]
    return and_join(listed_difference)


def show_diag_differences(repl: List[Tuple[int, int, int]], ins1: Dict[int, np.array], ins2: Dict[int, np.array], translation: Callable[[int], str] = str) -> str:
    repls_as_str = [(i, f"({seq_write_tuple[nuc1]})")
                    for i, nuc1, nuc2 in repl]
    ins_as_str = [
        (i, f"(insertion {''.join(seq_write_tuple[nuc] for nuc in frag)})") for i, frag in ins1.items()]
    del_as_str = [
        (i, f"(deletion {''.join(seq_write_tuple[nuc] for nuc in frag)})") for i, frag in ins2.items()]
    listed_difference = [f"{translation(i)} {desc}" for i, desc in sorted(
        repls_as_str + ins_as_str + del_as_str, key=lambda x: x[0])]
    return ", ".join(listed_difference)


def diag_textual_differences(repl: List[Tuple[int, int, int]], ins1: Dict[int, np.array], ins2: Dict[int, np.array], translation: Callable[[int], str] = str) -> str:
    repls_as_str = [(i, f"having a {seq_write_tuple[nuc1]}")
                    for i, nuc1, _ in repl]
    ins_as_str = [
        (i, f"having an insertion {''.join(seq_write_tuple[nuc] for nuc in frag)}") for i, frag in ins1.items()]
    del_as_str = [
        (i, f"having a deletion {''.join(seq_write_tuple[nuc] for nuc in frag)}") for i, frag in ins2.items()]
    listed_difference = [f"{desc} at position {translation(i)}" for i, desc in sorted(
        repls_as_str + ins_as_str + del_as_str, key=lambda x: x[0])]
    return and_join(listed_difference)


class DnaProcessor():

    def __init__(self, output_dir: str) -> None:
        self.infile: Optional[str] = None
        self.table: Optional[pd.DataFrame] = None
        self.aligned = False
        self.insertions = False
        self.relative_positions = False
        self.output_dir = output_dir

    def output(self, name) -> TextIO:
        filename = os.path.join(self.output_dir, name + ".txt")
        return open(filename, mode="w")

    def load_table(self, infile: str) -> None:
        with open(infile, errors='replace') as file:
            table = pd.read_csv(file, delimiter='\t').rename(
                columns=str.casefold).rename(columns=typos)
        if 'sequence' not in table.columns:
            raise ValueError("'sequences' or 'sequence' column is missing")
        if 'specimenid' not in table.columns:
            warnings.warn("Specimen IDs are not detected")
        else:
            table.set_index('specimenid', inplace=True)
        if len(table.columns) < 2:
            raise ValueError("'species' or another column need to be present")
        if self.aligned:
            table['sequence'] = table['sequence'].apply(Seq.from_str_notrim)
        else:
            table['sequence'] = table['sequence'].apply(Seq.from_str)
        self.table = table
        self.infile = infile

    def process_files(self, infile: str, reference_name: str, column: str, selection: List[str]) -> None:
        if not infile:
            raise ValueError('Input file is not given')
        if self.infile != infile:
            self.load_table(infile)
            column = "species"
            selection = []
        assert(self.table is not None)
        if len(selection) == 1:
            raise ValueError(
                "Please select at least two categories for comparison")
        reference_sequence = references[reference_name]
        alignment_displays = self.table['sequence'].apply(
            lambda seq: seq.align(reference_sequence)) if not self.aligned else self.table['sequence']
        if not self.insertions:
            self.table['sequence'].apply(Seq.reset_insertions)
        try:
            table = self.table.groupby(
                column)['sequence'].agg(combine_sequences)
        except ValueError as ex:
            raise ValueError("The sequences seem to not be aligned") from ex
        position_translator = table.iat[0].make_position_tranlator(
            reference_sequence)
        if selection:
            table = table.filter(items=selection).sort_index()
        with self.output("Aligments") as output:
            print("Alignments:", file=output)
            for specimen, alignment_str in alignment_displays.items():
                if not selection or self.table[column][specimen] in selection:
                    print(specimen, file=output)
                    print(alignment_str, file=output)
            output.write("\n")
        if self.relative_positions:
            self.report(table, column, reference_name,
                        lambda i: position_translator[i])
        else:
            self.report(table, column, reference_name)

    def report(self, table: pd.Series, column: str, reference_name: str, translation: Callable[[int], str] = str) -> None:
        with self.output("Difference_matrix") as matrixOutput, self.output("Difference_table") as tableOutput, self.output("Differences_description") as textOutput:
            print("", *table.index, sep='\t', file=matrixOutput)
            print(
                f"{column} 1\t{column} 2\treplacements\tinsertions 1\tinsertions 2", file=tableOutput)
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
                            repl, ins1, ins2, translation), sep='\t', file=tableOutput)
                    if difference_num > 0:
                        textFragments.append(
                            f"from {species2} in nucleotide {'position' if difference_num == 1 else 'positions'} " +
                            textual_differences(repl, ins1, ins2, translation))
                textOutput.write(
                    and_join(textFragments))
                textOutput.write("\n\n")
                matrixOutput.write("\n")

        diag_tableOutput = io.StringIO("")
        diag_textOutput = io.StringIO("")
        with self.output("Diagnostic_table") as diag_tableOutput, self.output("Diagnostics_description") as diag_textOutput:
            print(f"{column}\tUnique diagnostic differences",
                  file=diag_tableOutput)
            for species1 in table.index:
                other_species_seq: Seq = combine_sequences(table.drop(
                    labels=species1))
                repl, ins1, ins2 = differences(
                    table[species1], other_species_seq)
                text = show_diag_differences(repl, ins1, ins2, translation)
                if text:
                    print(species1, text, sep='\t', file=diag_tableOutput)
                else:
                    print(species1, "None", sep='\t', file=diag_tableOutput)
                text = diag_textual_differences(repl, ins1, ins2, translation)
                if text:
                    diag_textOutput.write(
                        f"{species1} differs from all other {column} in the dataset by ")
                    diag_textOutput.write(text)
                    diag_textOutput.write(
                        f" of the {reference_name} reference sequence.\n")
                else:
                    print(
                        f"{species1} has no unique diagnostic differences in comparison to the other {column} in the data set.", file=diag_textOutput)

    def choices(self) -> Dict[str, List[str]]:
        """
        Returns the columns that contain values for grouping
        """
        assert(self.table is not None)
        return {column_name: list(self.table[column_name].unique()) for column_name in self.table.columns if column_name != "sequence"}
