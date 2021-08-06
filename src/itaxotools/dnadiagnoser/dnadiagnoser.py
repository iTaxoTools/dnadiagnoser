#!/usr/bin/env python3

import sys
import tempfile


from .library.gui_utils import *
from .library.dnaprocessor import DnaProcessor
from .library.gui import DNADiagnoserGUI
from .library.resources import get_resource


def launch_gui() -> None:
    root = tk.Tk()

    def close_window():
        root.destroy()
        root.quit()

    root.title("DNAdiagnoser")
    if os.name == "nt":
        root.wm_iconbitmap(get_resource('dnadiagnoser.ico'))

    root.protocol("WM_DELETE_WINDOW", close_window)
    root.rowconfigure(0, weight=1)
    root.columnconfigure(0, weight=1)

    preview_dir = tempfile.mkdtemp()
    DNADiagnoserGUI(root, preview_dir=preview_dir)

    root.mainloop()
    root.quit()


def main() -> None:
    if len(sys.argv) >= 3:
        input = sys.argv[1]
        try:
            reference_name = sys.argv[3]
        except IndexError:
            reference_name = "Homo_sapiens_COI"
        output_dir = tempfile.mkdtemp()
        processor = DnaProcessor(output_dir)
        processor.process_files(input, reference_name, "species", [])
    else:
        launch_gui()


if __name__ == "__main__":
    main()
