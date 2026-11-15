import argparse
import csv
import sys

# Código genético estándar (ADN -> Aminoácido de 1 letra)
CODON_TABLE = {
    "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K",
    "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
    "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",
    "TAC": "Y", "TAT": "Y", "TAA": "*", "TAG": "*",
    "TGC": "C", "TGT": "C", "TGA": "*", "TGG": "W",
}

# Clasificación bioquímica de cadenas laterales (grupos R)
AA_PROPERTIES = {
    "A": "Nonpolar", "V": "Nonpolar", "L": "Nonpolar", "I": "Nonpolar",
    "P": "Nonpolar", "F": "Nonpolar", "W": "Nonpolar", "M": "Nonpolar",
    "G": "Nonpolar",
    "S": "Polar", "T": "Polar", "C": "Polar", "Y": "Polar",
    "N": "Polar", "Q": "Polar",
    "K": "Positively Charged", "R": "Positively Charged", "H": "Positively Charged",
    "D": "Negatively Charged", "E": "Negatively Charged",
    "*": "Stop Codon"
}


def parse_fasta(filepath):
    """Lee un archivo FASTA y devuelve una lista de tuplas (header, sequence)."""
    sequences = []
    header = None
    seq_chunks = []

    try:
        with open(filepath, "r", encoding="utf-8") as file:
            for line in file:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if header:
                        sequences.append((header, "".join(seq_chunks).upper()))
                    header = line[1:]
                    seq_chunks = []
                else:
                    seq_chunks.append(line)
            if header:
                sequences.append((header, "".join(seq_chunks).upper()))
    except FileNotFoundError:
        print(f"[ERROR] No se encontró el archivo: {filepath}")
        sys.exit(1)

    return sequences


def translate_dna(dna_seq):
    """Traduce una secuencia de ADN (múltiplo de 3) a secuencia proteica."""
    protein = []
    for i in range(0, len(dna_seq) - (len(dna_seq) % 3), 3):
        codon = dna_seq[i:i+3]
        protein.append(CODON_TABLE.get(codon, "X"))
    return "".join(protein)


def classify_variant(ref_codon, sample_codon):
    """Clasifica una mutación puntual entre dos codones y devuelve sus metadatos."""
    ref_aa = CODON_TABLE.get(ref_codon, "X")
    sample_aa = CODON_TABLE.get(sample_codon, "X")

    if ref_aa == sample_aa:
        return ref_aa, sample_aa, "SYNONYMOUS", "No amino acid change"
    if sample_aa == "*":
        return ref_aa, sample_aa, "NONSENSE", "Premature STOP codon introduced"

    ref_prop = AA_PROPERTIES.get(ref_aa, "Unknown")
    sample_prop = AA_PROPERTIES.get(sample_aa, "Unknown")
    if ref_prop == sample_prop:
        impact = f"Conservative ({ref_prop} -> {sample_prop})"
    else:
        impact = f"Non-conservative ({ref_prop} -> {sample_prop})"

    return ref_aa, sample_aa, "MISSENSE", impact


def analyze_sample(ref_name, ref_seq, sample_name, sample_seq):
    """Compara una muestra contra la referencia y devuelve una lista de diccionarios con variantes."""
    if len(ref_seq) != len(sample_seq):
        print(f"[WARNING] Saltando '{sample_name}': longitud distinta ({len(sample_seq)} bp vs {len(ref_seq)} bp).")
        return []

    ref_prot = translate_dna(ref_seq)
    sample_prot = translate_dna(sample_seq)

    print("=" * 70)
    print(f"SAMPLE: {sample_name} (vs {ref_name})")
    print(f"Ref Protein   : {ref_prot}")
    print(f"Sample Protein: {sample_prot}")
    print("-" * 70)

    detected_variants = []

    for codon_idx in range(0, len(ref_seq), 3):
        ref_codon = ref_seq[codon_idx:codon_idx+3]
        sample_codon = sample_seq[codon_idx:codon_idx+3]

        if ref_codon != sample_codon:
            aa_pos = (codon_idx // 3) + 1
            ref_aa, sample_aa, v_type, impact = classify_variant(ref_codon, sample_codon)

            variant_record = {
                "reference": ref_name,
                "sample": sample_name,
                "aa_position": aa_pos,
                "ref_codon": ref_codon,
                "sample_codon": sample_codon,
                "ref_aa": ref_aa,
                "sample_aa": sample_aa,
                "mutation_type": v_type,
                "biochemical_impact": impact
            }
            detected_variants.append(variant_record)

            print(f"  Posición AA {aa_pos:>2} | {ref_codon} -> {sample_codon} | AA: {ref_aa} -> {sample_aa} | [{v_type}] {impact}")

    if not detected_variants:
        print("  Sin variantes detectadas.")
    print()

    return detected_variants


def export_to_csv(variants, output_path):
    """Exporta la lista de variantes detectadas a un archivo CSV."""
    fieldnames = [
        "reference", "sample", "aa_position", "ref_codon",
        "sample_codon", "ref_aa", "sample_aa", "mutation_type", "biochemical_impact"
    ]
    with open(output_path, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(variants)
    print(f"[OK] Reporte exportado exitosamente a: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Traduce secuencias codificantes (FASTA) y clasifica mutaciones puntuales."
    )
    parser.add_argument(
        "-i", "--input",
        default="data/sample_sequences.fasta",
        help="Ruta al archivo FASTA de entrada (por defecto: data/sample_sequences.fasta)"
    )
    parser.add_argument(
        "-o", "--output",
        help="Ruta opcional para exportar los resultados en formato CSV (ej. results.csv)"
    )
    args = parser.parse_args()

    records = parse_fasta(args.input)
    if len(records) < 2:
        print("[ERROR] El archivo FASTA debe contener al menos 2 secuencias (Referencia y al menos 1 Muestra).")
        sys.exit(1)

    ref_header, ref_dna = records[0]
    all_variants = []

    for sample_header, sample_dna in records[1:]:
        sample_variants = analyze_sample(ref_header, ref_dna, sample_header, sample_dna)
        all_variants.extend(sample_variants)

    if args.output:
        export_to_csv(all_variants, args.output)


if __name__ == "__main__":
    main()
