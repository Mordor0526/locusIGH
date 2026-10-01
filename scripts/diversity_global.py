#!/usr/bin/env python3
"""
diversity_global.py
Calculate nucleotide diversity indices per species from an IgH locus MSA.

Computes: n, S (segregating sites), h (haplotypes), Hd (haplotype diversity),
pi (nucleotide diversity), Tajima's D.

Usage:
    python3 diversity_global.py IgH_103_aligned.msa --outgroup Camelus_bactrianus \
        -o diversity_results.csv

Species assignment is based on the last character of each sequence ID:
  g = L. guanicoe, l = L. glama, v = V. vicugna, a = V. pacos
"""

import argparse
import sys
import numpy as np
from Bio import SeqIO

SPECIES_MAP = {
    "g": "L_guanicoe",
    "l": "L_glama",
    "v": "V_vicugna",
    "a": "V_pacos",
}


def calc_diversity(subset):
    """Calculate diversity indices for a list of (name, seq) tuples."""
    names, seqs = zip(*subset)
    n = len(seqs)
    L = len(seqs[0])

    # Convert to numpy arrays
    arr = np.array([list(s) for s in seqs], dtype="U1")

    # Numeric encoding: A=0, C=1, G=2, T=3, other=99
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

    # Haplotypes
    seq_list = list(seqs)
    haps = set(seq_list)
    h = len(haps)
    Hd = (n / (n - 1)) * (
        1 - sum((seq_list.count(hp) / n) ** 2 for hp in haps)
    ) if n > 1 else 0

    # Nucleotide diversity (pi) - vectorized pairwise
    total_diffs = 0
    total_valid = 0
    for i in range(n):
        for j in range(i + 1, n):
            both_valid = (num[i] < 99) & (num[j] < 99)
            total_valid += np.sum(both_valid)
            total_diffs += np.sum(both_valid & (num[i] != num[j]))
    pi = total_diffs / total_valid if total_valid > 0 else 0

    # Tajima's D
    a1 = sum(1.0 / i for i in range(1, n))
    a2 = sum(1.0 / (i ** 2) for i in range(1, n))
    b1 = (n + 1) / (3 * (n - 1))
    b2 = 2 * (n ** 2 + n + 3) / (9 * n * (n - 1))
    c1 = b1 - 1 / a1
    c2 = b2 - (n + 2) / (a1 * n) + a2 / (a1 ** 2)
    e1 = c1 / a1
    e2 = c2 / (a1 ** 2 + a2)
    k_hat = 2 * total_diffs / (n * (n - 1))
    tajD_num = k_hat - S / a1
    tajD_den = (e1 * S + e2 * S * (S - 1)) ** 0.5 if S > 0 else 1
    tajD = tajD_num / tajD_den if tajD_den > 0 else 0

    return n, S, h, Hd, pi, tajD


def main():
    parser = argparse.ArgumentParser(
        description="Nucleotide diversity per species for IgH locus"
    )
    parser.add_argument("input", help="Input MSA file (FASTA)")
    parser.add_argument("--outgroup", default="Camelus_bactrianus",
                        help="Outgroup sequence ID to exclude (default: Camelus_bactrianus)")
    parser.add_argument("-o", "--output", default=None,
                        help="Output CSV file (default: stdout)")
    args = parser.parse_args()

    # Load sequences
    recs = [
        (r.id, str(r.seq).upper())
        for r in SeqIO.parse(args.input, "fasta")
        if r.id != args.outgroup
    ]
    print(f"Loaded {len(recs)} ingroup sequences", file=sys.stderr)

    # Group by species
    groups = {}
    for name, seq in recs:
        sp = SPECIES_MAP.get(name[-1], "unknown")
        groups.setdefault(sp, []).append((name, seq))

    # Output
    out = open(args.output, "w") if args.output else sys.stdout
    out.write("group,n,S,h,Hd,pi,TajD\n")

    for label, subset in [("ALL_CSA", recs)] + [
        (sp, groups[sp])
        for sp in ["L_guanicoe", "L_glama", "V_vicugna", "V_pacos"]
        if sp in groups
    ]:
        n, S, h, Hd, pi, tajD = calc_diversity(subset)
        out.write(f"{label},{n},{S},{h},{Hd:.5f},{pi:.5f},{tajD:.4f}\n")
        out.flush()
        print(f"  {label} done (n={n})", file=sys.stderr)

    if args.output:
        out.close()
        print(f"Results written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
