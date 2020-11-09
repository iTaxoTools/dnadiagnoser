import os.path
import tkinter as tk
import tkinter.filedialog
import tkinter.ttk as ttk
from typing import Literal, Any, Dict, Tuple, List


class FileChooser():
    """
    Creates a frame with a label, entry and browse button for choosing files
    """

    def __init__(self, parent: Any, *, label: str, mode: Literal["open", "save"]):
        self.frame = ttk.Frame(parent)
        self.frame.columnconfigure([0, 1], weight=1)
        self.label = ttk.Label(self.frame, text=label)
        self.file_var = tk.StringVar()
        self.entry = ttk.Entry(self.frame, textvariable=self.file_var)
        if mode == "open":
            self._dialog = tk.filedialog.askopenfilename
        elif mode == "save":
            self._dialog = tk.filedialog.asksaveasfilename

        def browse() -> None:
            if (newpath := self._dialog()):
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

    def __init__(self, parent: tk.Misc, *, height: int, selectmode: Literal["browse", "extended"], values: List[str]) -> None:
        self.list = values
        self.listbox = tk.Listbox(
            parent, height=height, selectmode=selectmode, listvariable=tk.StringVar(value=values))
        self.grid = self.listbox.grid

    def selection(self) -> List[str]:
        """
        Returns the list of items that are currently selected
        """
        return [self.list[i] for i in self.listbox.curselection()]


class ColumnSelector():

    def __init__(self, parent: tk.Misc, columns: Dict[str, List[str]]) -> None:
        self.notebook = ttk.Notebook(parent)
        self.lists = []
        for column_name in columns:
            self.lists.append(Listbox(self.notebook, height=10,
                                      values=columns[column_name], selectmode="extended"))
            self.notebook.add(self.lists[-1].listbox, text=column_name)
        self.grid = self.notebook.grid

    def selection(self) -> Tuple[str, List[str]]:
        i = self.notebook.index("current")
        return (self.notebook.tab(i)["text"], self.lists[i].selection())
