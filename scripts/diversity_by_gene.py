#!/usr/bin/env python3
"""
diversity_by_gene.py
Calculate nucleotide diversity (pi) and segregating sites (S) per annotated
gene region for each species in the IgH locus MSA.

Requires:
  - Aligned MSA (FASTA)
  - Gene coordinates file (TSV with columns: gene, aln_start, aln_end, type)

Usage:
    python3 diversity_by_gene.py IgH_103_aligned.msa gene_coords_aln.tsv \
        --outgroup Camelus_bactrianus -o diversity_by_gene.csv
"""

import argparse
import csv
import sys
import numpy as np
from Bio import SeqIO

SPECIES_MAP = {
    "g": "L_guanicoe",
    "l": "L_glama",
    "v": "V_vicugna",
    "a": "V_pacos",
}


def diversity_region(subset, start, end):
    """Calculate n, length, S, and pi for a specific alignment region."""
    names, seqs = zip(*subset)
    n = len(seqs)
    L = end - start
    region = [s[start:end] for s in seqs]

    arr = np.array([list(r) for r in region], dtype="U1")
    mapping = {"A": 0, "C": 1, "G": 2, "T": 3}
    num = np.full((n, L), 99, dtype=np.int8)
    for base, val in mapping.items():
        num[arr == base] = val

    # Segregating sites
    S = 0
    for j in range(L):
        col = set(arr[:, j]) - {"N", "-"}
        if len(col) > 1:
            S += 1

    # Pi - pairwise
    total_diffs = 0
    total_valid = 0
    for i in range(n):
        for j in range(i + 1, n):
            bv = (num[i] < 99) & (num[j] < 99)
            total_valid += np.sum(bv)
            total_diffs += np.sum(bv & (num[i] != num[j]))
    pi = total_diffs / total_valid if total_valid > 0 else 0

    return n, L, S, pi


def main():
    parser = argparse.ArgumentParser(
        description="Diversity per gene region in the IgH locus"
    )
    parser.add_argument("input", help="Input MSA file (FASTA)")
    parser.add_argument("coords", help="Gene coordinates TSV (gene, aln_start, aln_end, type)")
    parser.add_argument("--outgroup", default="Camelus_bactrianus",
                        help="Outgroup ID to exclude")
    parser.add_argument("-o", "--output", default=None, help="Output CSV (default: stdout)")
    args = parser.parse_args()

    # Load sequences
    recs = [
        (r.id, str(r.seq).upper())
        for r in SeqIO.parse(args.input, "fasta")
        if r.id != args.outgroup
    ]
    print(f"Loaded {len(recs)} ingroup sequences", file=sys.stderr)

    # Load gene coords
    genes = []
    with open(args.coords) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            genes.append((row["gene"], int(row["aln_start"]), int(row["aln_end"]), row["type"]))
    print(f"Loaded {len(genes)} gene regions", file=sys.stderr)

    # Group by species
    groups = {"ALL_CSA": recs}
    for name, seq in recs:
        sp = SPECIES_MAP.get(name[-1], "unknown")
        groups.setdefault(sp, []).append((name, seq))

    # Output
    out = open(args.output, "w") if args.output else sys.stdout
    out.write("gene,type,length,group,n,S,pi\n")

    for gname, start, end, gtype in genes:
        for grp in ["ALL_CSA", "L_guanicoe", "L_glama", "V_vicugna", "V_pacos"]:
            if grp in groups:
                n, L, S, pi = diversity_region(groups[grp], start - 1, end)
                out.write(f"{gname},{gtype},{L},{grp},{n},{S},{pi:.6f}\n")
                out.flush()
        print(f"  {gname} done", file=sys.stderr)

    if args.output:
        out.close()
        print(f"Results written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
