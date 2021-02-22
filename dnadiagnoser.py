#!/usr/bin/env python3

import sys
import tempfile


from library.gui_utils import *
from library.dnaprocessor import DnaProcessor
from library.gui import DNADiagnoserGUI


def launch_gui() -> None:
    root = tk.Tk()

    def close_window():
        root.destroy()
        root.quit()

    root.title("DNAdiagnoser")
    if os.name == "nt":
        root.wm_iconbitmap(os.path.join(
            sys.path[0], 'morphometricanalyzer.ico'))

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
