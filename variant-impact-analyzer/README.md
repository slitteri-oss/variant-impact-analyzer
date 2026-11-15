# Variant Impact Analyzer

A lightweight Python command-line tool for translating coding DNA sequences (CDS) and classifying point mutations between a reference sequence and a sample sequence in FASTA format.

## Features
- DNA to Protein Translation: Translates aligned coding sequences using the standard genetic code.
- Variant Classification: Identifies codon changes and classifies mutations into Synonymous (Silent), Missense (evaluating side-chain chemical conservation), and Nonsense (premature stop codon).

## Project Structure
- data/sample_sequences.fasta: Example reference and mutated CDS
- src/analyzer.py: Core translation and variant classification logic
- .gitignore: Ignored cache and local files
- README.md: Project documentation

## Usage
Run the analyzer on the included sample dataset with: python3 src/analyzer.py
