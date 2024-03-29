import os.path
import tkinter as tk
import tkinter.filedialog as tkfiledialog
import tkinter.ttk as ttk
from typing import Any, Dict, Tuple, List, Optional


class FileChooser():
    """
    Creates a frame with a label, entry and browse button for choosing files
    """

    def __init__(self, parent: Any, *, label: str, mode: str):
        self.frame = ttk.Frame(parent)
        self.frame.columnconfigure([0, 1], weight=1)
        self.label = ttk.Label(self.frame, text=label)
        self.file_var = tk.StringVar()
        self.entry = ttk.Entry(self.frame, textvariable=self.file_var)
        if mode == "open":
            self._dialog = tkfiledialog.askopenfilename
        elif mode == "save":
            self._dialog = tkfiledialog.asksaveasfilename

        def browse() -> None:
            newpath: Optional[str] = self._dialog()
            if newpath:
                try:
                    newpath = os.path.relpath(newpath)
                except:
                    newpath = os.path.abspath(newpath)
                self.file_var.set(newpath)

        self.button = ttk.Button(self.frame, text="Browse", command=browse)

        self.label.grid(row=0, column=0, sticky='nws')
        self.entry.grid(row=1, column=0, sticky='nwse')
        self.button.grid(row=1, column=1)
        self.grid = self.frame.grid


class LabeledEntry():
    """
    Group of a label, entry and a string variable
    """

    def __init__(self, parent: tk.Misc, *, label: str):
        self.frame = ttk.Frame(parent)
        self.label = ttk.Label(self.frame, text=label)
        self.var = tk.StringVar()
        self.entry = ttk.Entry(self.frame, textvariable=self.var)
        self.frame.columnconfigure(1, weight=1)
        self.label.grid(column=0, row=0)
        self.entry.grid(column=1, row=0)
        self.grid = self.frame.grid


class LabeledCombobox():
    """
    Group of a label, Combobox and a string variable
    """

    def __init__(self, parent: tk.Misc, *, label: str, values: List[str], readonly: bool):
        self.frame = ttk.Frame(parent)
        self.label = ttk.Label(self.frame, text=label)
        self.var = tk.StringVar()
        self.combobox = ttk.Combobox(
            self.frame, textvariable=self.var, values=values)
        if readonly:
            self.combobox.configure(state='readonly')
            if values:
                self.combobox.current(0)
        self.frame.columnconfigure(1, weight=1)
        self.label.grid(column=0, row=0)
        self.combobox.grid(column=1, row=0)
        self.grid = self.frame.grid


class Listbox():
    """
    Wrapper for a read-only tk.Listbox with a method that returns the selection
    """

    def __init__(self, parent: tk.Misc, *, height: int, selectmode: str, values: List[str]) -> None:
        self.list = values
        listvar = tk.StringVar(value=" ".join(value.replace(' ', r'\ ') for value in values))
        self.listbox = tk.Listbox(
            parent, height=height, selectmode=selectmode, listvariable=listvar)
        self.grid = self.listbox.grid

    def selection(self) -> List[str]:
        """
        Returns the list of items that are currently selected
        """
        return [self.list[i] for i in self.listbox.curselection()]


class ColumnSelector():

    def __init__(self, parent: tk.Misc) -> None:
        self.notebook = ttk.Notebook(parent)
        self.lists: List[Listbox] = []
        self.grid = self.notebook.grid

    def set_columns(self, columns: Dict[str, List[str]]) -> None:
        for i in range(self.notebook.index("end")):
            self.notebook.forget(i)
        self.lists = []
        for column_name in columns:
            self.lists.append(Listbox(self.notebook, height=10,
                                      values=columns[column_name], selectmode="extended"))
            self.notebook.add(self.lists[-1].listbox, text=column_name)

    def selection(self) -> Optional[Tuple[str, List[str]]]:
        if self.notebook.index("end") == 0:
            return None
        i = self.notebook.index("current")
        return (self.notebook.tab(i)["text"], self.lists[i].selection())
