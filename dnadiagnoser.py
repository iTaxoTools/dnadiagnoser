#!/usr/bin/env python3

import sys
import tkinter.messagebox
import warnings
from typing import Optional


from library.gui_utils import *
from library.dnaprocessor import DnaProcessor, references


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

    options_frame = ttk.Frame(root)
    aligned_var = tk.BooleanVar(value=False)
    relative_positions_var = tk.BooleanVar(value=False)

    def set_relative_positions() -> None:
        processor.relative_positions = relative_positions_var.get()
    relative_positions_btn = ttk.Checkbutton(
        options_frame, variable=relative_positions_var, offvalue=False, onvalue=True, text="Print positions relative to the reference sequence", command=set_relative_positions)
    relative_positions_btn.state(['disabled'])

    def set_aligned() -> None:
        processor.aligned = aligned_var.get()
        processor.relative_positions &= processor.aligned
        relative_positions_btn.state([
            '!disabled' if processor.aligned else 'disabled'])
    aligned_btn = ttk.Checkbutton(
        options_frame, variable=aligned_var, offvalue=False, onvalue=True, text="Already aligned", command=set_aligned)
    aligned_btn.grid(row=0, column=0, sticky='w')
    relative_positions_btn.grid(row=1, column=0, sticky='w')

    def load() -> None:
        try:
            infile: Optional[str] = inputchooser.file_var.get()
            if infile:
                processor.load_table(infile)
                column_selector.set_columns(processor.choices())
        except Exception as ex:
            tkinter.messagebox.showerror("Error", str(ex))

    def process() -> None:
        try:
            with warnings.catch_warnings(record=True) as warns:
                selection = column_selector.selection()
                if selection:
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

    options_frame.grid(row=2, column=2)

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
