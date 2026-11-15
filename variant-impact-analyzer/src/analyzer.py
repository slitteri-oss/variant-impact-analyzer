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

# Clasificación bioquímica básica de cadenas laterales (grupos R)
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

    return sequences


def translate_dna(dna_seq):
    """Traduce una secuencia de ADN (múltiplo de 3) a secuencia proteica."""
    protein = []
    for i in range(0, len(dna_seq) - (len(dna_seq) % 3), 3):
        codon = dna_seq[i:i+3]
        protein.append(CODON_TABLE.get(codon, "X"))
    return "".join(protein)


def analyze_variants(ref_name, ref_seq, sample_name, sample_seq):
    """Compara referencia vs muestra codón por codón e imprime el reporte."""
    if len(ref_seq) != len(sample_seq):
        print(f"[ERROR] Longitudes distintas: Ref ({len(ref_seq)} bp) vs Sample ({len(sample_seq)} bp).")
        print("Este analizador MVP requiere secuencias codificantes alineadas de igual longitud.")
        sys.exit(1)

    ref_prot = translate_dna(ref_seq)
    sample_prot = translate_dna(sample_seq)

    print("=" * 65)
    print(f"VARIANT IMPACT REPORT: {sample_name} vs {ref_name}")
    print("=" * 65)
    print(f"Ref Protein   : {ref_prot}")
    print(f"Sample Protein: {sample_prot}")
    print("-" * 65)

    variants_found = 0

    for codon_idx in range(0, len(ref_seq), 3):
        ref_codon = ref_seq[codon_idx:codon_idx+3]
        sample_codon = sample_seq[codon_idx:codon_idx+3]

        if ref_codon != sample_codon:
            variants_found += 1
            aa_pos = (codon_idx // 3) + 1
            ref_aa = CODON_TABLE.get(ref_codon, "X")
            sample_aa = CODON_TABLE.get(sample_codon, "X")

            if ref_aa == sample_aa:
                v_type = "SYNONYMOUS (Silent)"
                impact = "No amino acid change"
            elif sample_aa == "*":
                v_type = "NONSENSE (Truncating)"
                impact = "Premature STOP codon introduced!"
            else:
                v_type = "MISSENSE"
                ref_prop = AA_PROPERTIES.get(ref_aa, "Unknown")
                sample_prop = AA_PROPERTIES.get(sample_aa, "Unknown")
                if ref_prop == sample_prop:
                    impact = f"Conservative ({ref_prop} -> {sample_prop})"
                else:
                    impact = f"Non-conservative ({ref_prop} -> {sample_prop})"

            print(f"Posición AA {aa_pos:>2} | Codón: {ref_codon} -> {sample_codon} | AA: {ref_aa} -> {sample_aa}")
            print(f"   Tipo   : {v_type}")
            print(f"   Impacto: {impact}\n")

    if variants_found == 0:
        print("No se detectaron variantes entre las secuencias.")
    print("=" * 65)


if __name__ == "__main__":
    fasta_path = "data/sample_sequences.fasta"
    records = parse_fasta(fasta_path)

    if len(records) < 2:
        print("[ERROR] El archivo FASTA debe contener al menos 2 secuencias (Referencia y Muestra).")
        sys.exit(1)

    ref_header, ref_dna = records[0]
    sample_header, sample_dna = records[1]

    analyze_variants(ref_header, ref_dna, sample_header, sample_dna)
