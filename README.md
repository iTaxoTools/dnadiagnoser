# dnadiagnoser
Identifies diagnostic DNA nucleotides

## Description
This program generates diagnostic information from a tab-file of sequences

### Input file
The input should be a tab-file with DNA sequences.
A `specimenid` (or a variation of this) column is expected, but not required.
A `sequences` column (or a variation) is required.
At least one other column is required.

## Usage

The interface contains two main buttons.

The `Load` button allows to load the input file and select the column and values that will be processed.

The `Process` button processes the table and writes the output.
If the current input file has not been loaded, the column `species` is selected and all the values are processed.
If the file has been loaded, but not values have been selected, then all the values of the selected column are processed.
Otherwise, only the selected values are processed.

The `Already aligned` checkbox can be used, if the sequences in the input file are already aligned and have the same length.
Then the step of aligning them with the reference sequence is skipped and the sequences are compared as they are in the input file.

The `Print position relative to the reference sequence` checkbox enables relative position translation.
The first sequence in the input file is aligned with the chosen reference sequence to calculate the correspondence between positions in input sequences and the reference sequence.
Then this correspondence is used to translate the position written to the output file.
The insertions relative to the reference sequence are written as `n+i`, which means `i`\-th nucleotide inserted after `n`\-th nucleotide of the reference sequence.

### Command-line mode

Dnadiagnoser include a rudimentary command-line mode. It can be used in the following way:
```
dnadiagnoser input_file output_directory [reference_sequence_name]
```

It will process the `input_file` (which should be in the Tab format)
and write the output into the directory specified by the path `output_directory`.
`reference_sequence_name` is an optional argument that specifies the name of the reference sequence.

## Reference sequences
The file `src/itaxotools/dnadiagnoser/resources/reference_sequences.tab` contains the reference sequences used for alignment.
Each lines has the format:
```
sequence_name<Tab>sequence
```

## Scores for alignment
The file `src/itaxotools/dnadiagnoser/resources/scores.tab` contains the scores used in the sequence alignment.
Each line has the format:
```
score_identifier<Tab>value
```

The possible scores are:
* `gap penalty`: Score to open a gap in the middle of a sequence
* `gap extend penalty`: Score to extend an existing gap in the middle of a sequence
* `end gap penalty`: Score to create a gap at an end of a sequence.
* `end gap extend penalty`: Score to extend a gap at an end of a sequence.
* `match score`: Score for matching nucleotides
* `mismatch score`: Score for non-matching nucleotides
