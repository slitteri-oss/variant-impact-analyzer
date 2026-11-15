# Variant Impact Analyzer

A Python command-line tool for translating coding DNA sequences (CDS), comparing multiple samples against a reference sequence in FASTA format, and classifying point mutations along with their biochemical impact.

## Features
- Multi-sample FASTA comparison against a reference coding sequence.
- DNA to protein translation using the standard genetic code.
- Point mutation classification: Synonymous (silent), Missense (evaluating side-chain chemical conservation), and Nonsense (premature stop codon).
- Command-line interface (CLI) built with argparse and structured CSV report export.
- Unit test suite covering translation and variant classification logic.

## Project Structure
- data/sample_sequences.fasta: Example reference and mutated CDS dataset
- src/analyzer.py: Core CLI, translation, variant classification, and CSV export logic
- tests/test_analyzer.py: Unit tests for sequence translation and mutation classification
- .gitignore: Ignored cache, temporary CSV outputs, and local files
- README.md: Project documentation

## Usage
Run analysis on the default dataset:
python3 src/analyzer.py

Specify a custom input FASTA and export results to CSV:
python3 src/analyzer.py -i data/sample_sequences.fasta -o variants_report.csv

Run the unit test suite:
python3 -m unittest discover -s tests -v
