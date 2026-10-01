#!/usr/bin/env python3
"""
Tajima's D by gene region of the IGH locus in South American camelids.
Usage: python3 tajima_d_igh.py
Expects: IgH_103_aligned.msa and gene_coords_aln.tsv in the working directory.
"""

import sys, math, os
from collections import defaultdict

# ── Configuration ──────────────────────────────────────────────────
ALN_FILE = os.path.expanduser("~/locusigh/IgH_103_aligned.msa")
COORDS_FILE = os.path.expanduser("~/locusigh/gene_coords_aln.tsv")

# Species assignment (same logic as previous scripts)
SPECIES_MAP = {}  # filled when reading alignment

def assign_species(name):
    n = name.lower()
    if 'cbac' in n or 'bactrianus' in n or n.startswith('cb'):
        return 'C_bactrianus'
    # guanaco ids: numeric + 'g'
    if n.endswith('g') and n[:-1].isdigit():
        return 'L_guanicoe'
    if n.endswith('l') and n[:-1].isdigit():
        return 'L_glama'
    if n.endswith('v') and n[:-1].isdigit():
        return 'V_vicugna'
    if n.endswith('a') and n[:-1].isdigit():
        return 'V_pacos'
    # fallback patterns
    for tag, sp in [('guan', 'L_guanicoe'), ('glam', 'L_glama'), ('llam', 'L_glama'),
                    ('vicu', 'V_vicugna'), ('paco', 'V_pacos'), ('alpa', 'V_pacos')]:
        if tag in n:
            return sp
    return 'unknown'


def read_fasta(path):
    """Read FASTA alignment, return list of (name, sequence)."""
    seqs = []
    name, buf = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if name is not None:
                    seqs.append((name, ''.join(buf)))
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.upper())
    if name is not None:
        seqs.append((name, ''.join(buf)))
    return seqs


def read_gene_coords(path):
    """Read gene_coords_aln.tsv → list of (gene, start, end, type)."""
    regions = []
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                regions.append((parts[0], int(parts[1]), int(parts[2]), parts[3]))
    return regions


def extract_region(seqs, start, end):
    """Extract sub-alignment for a region [start, end] (0-based inclusive)."""
    return [(name, seq[start:end+1]) for name, seq in seqs]


def count_segregating_sites(sub_seqs):
    """Count segregating sites (S) in a sub-alignment. Returns S and n (sample size)."""
    n = len(sub_seqs)
    if n < 4:
        return 0, n, 0  # need at least 4 for Tajima's D
    seqlen = len(sub_seqs[0][1])
    S = 0
    usable_sites = 0
    for j in range(seqlen):
        bases = []
        for i in range(n):
            b = sub_seqs[i][1][j]
            if b in 'ACGT':
                bases.append(b)
        if len(bases) < n * 0.5:  # skip sites with >50% missing
            continue
        usable_sites += 1
        if len(set(bases)) > 1:
            S += 1
    return S, n, usable_sites


def tajima_d(sub_seqs):
    """
    Calculate Tajima's D for a set of aligned sequences.
    Returns (D, S, pi, n, n_sites) or None if insufficient data.
    """
    n = len(sub_seqs)
    if n < 4:
        return None

    seqlen = len(sub_seqs[0][1])

    # Calculate pi (average pairwise differences) and S simultaneously
    S = 0
    total_diffs = 0
    n_pairs = n * (n - 1) / 2
    usable_sites = 0

    for j in range(seqlen):
        bases = []
        for i in range(n):
            b = sub_seqs[i][1][j]
            if b in 'ACGT':
                bases.append(b)
        nn = len(bases)
        if nn < n * 0.5:
            continue
        usable_sites += 1
        alleles = set(bases)
        if len(alleles) > 1:
            S += 1
            # Count pairwise differences at this site
            counts = defaultdict(int)
            for b in bases:
                counts[b] += 1
            site_diffs = 0
            allele_list = list(counts.keys())
            for a in range(len(allele_list)):
                for b in range(a+1, len(allele_list)):
                    site_diffs += counts[allele_list[a]] * counts[allele_list[b]]
            # Normalise by actual pairs at this site
            site_pairs = nn * (nn - 1) / 2
            total_diffs += site_diffs / site_pairs

    if S == 0 or usable_sites == 0:
        return {'D': float('nan'), 'S': S, 'pi': 0, 'theta_W': 0,
                'n': n, 'n_sites': usable_sites, 'pi_per_site': 0}

    pi = total_diffs  # total pairwise differences (normalised per site per pair, summed across sites)

    # Watterson's theta
    a1 = sum(1.0/i for i in range(1, n))
    a2 = sum(1.0/(i*i) for i in range(1, n))
    theta_W = S / a1

    # Tajima's D coefficients
    b1 = (n + 1) / (3 * (n - 1))
    b2 = 2 * (n*n + n + 3) / (9 * n * (n - 1))
    c1 = b1 - 1/a1
    c2 = b2 - (n + 2) / (a1 * n) + a2 / (a1 * a1)
    e1 = c1 / a1
    e2 = c2 / (a1*a1 + a2)

    d = pi - theta_W
    var_d = e1 * S + e2 * S * (S - 1)

    if var_d <= 0:
        D = float('nan')
    else:
        D = d / math.sqrt(var_d)

    pi_per_site = pi / usable_sites if usable_sites > 0 else 0

    return {'D': D, 'S': S, 'pi': pi, 'theta_W': theta_W,
            'n': n, 'n_sites': usable_sites, 'pi_per_site': pi_per_site}


def main():
    print("Reading alignment...")
    seqs = read_fasta(ALN_FILE)
    print(f"  {len(seqs)} sequences, {len(seqs[0][1])} sites")

    # Assign species (exclude outgroup for within-SAC analyses)
    species_groups = defaultdict(list)
    sac_seqs = []
    for name, seq in seqs:
        sp = assign_species(name)
        if sp == 'C_bactrianus':
            continue
        species_groups[sp].append((name, seq))
        sac_seqs.append((name, seq))

    print(f"\nSAC samples: {len(sac_seqs)}")
    for sp in sorted(species_groups):
        print(f"  {sp}: n={len(species_groups[sp])}")

    # Read gene coordinates
    regions = read_gene_coords(COORDS_FILE)
    print(f"\nGene regions: {len(regions)}")

    # Define broader functional categories
    func_categories = {
        'VHH': [],
        'VH_conventional': [],
        'D_segments': [],
        'J_segments': [],
        'C_HCAb': [],
        'C_conventional': [],
        'C_other': [],
    }
    for gene, start, end, gtype in regions:
        if 'vhh' in gene.lower():
            func_categories['VHH'].append((gene, start, end))
        elif gtype == 'VH':
            func_categories['VH_conventional'].append((gene, start, end))
        elif gtype == 'D':
            func_categories['D_segments'].append((gene, start, end))
        elif gtype == 'J':
            func_categories['J_segments'].append((gene, start, end))
        elif gtype == 'C_HCAb':
            func_categories['C_HCAb'].append((gene, start, end))
        elif gtype == 'C_conv':
            func_categories['C_conventional'].append((gene, start, end))
        elif gtype == 'C':
            func_categories['C_other'].append((gene, start, end))

    # ── Analysis 1: Tajima's D per individual gene, all SAC ──
    print("\n" + "="*80)
    print("TAJIMA'S D PER GENE REGION — ALL SAC (n={})".format(len(sac_seqs)))
    print("="*80)
    print(f"\n{'Gene':<16} {'Type':<10} {'Length':>7} {'S':>6} {'pi':>10} {'theta_W':>10} {'D':>8}")
    print("-"*80)

    gene_results = []
    for gene, start, end, gtype in regions:
        sub = extract_region(sac_seqs, start, end)
        result = tajima_d(sub)
        if result:
            print(f"{gene:<16} {gtype:<10} {end-start+1:>7} {result['S']:>6} "
                  f"{result['pi']:>10.4f} {result['theta_W']:>10.4f} {result['D']:>8.4f}")
            gene_results.append((gene, gtype, end-start+1, result))

    # ── Analysis 2: Tajima's D per functional category, all SAC ──
    print("\n" + "="*80)
    print("TAJIMA'S D PER FUNCTIONAL CATEGORY — ALL SAC")
    print("="*80)
    print(f"\n{'Category':<20} {'Genes':>5} {'Length':>7} {'S':>6} {'pi':>10} {'theta_W':>10} {'D':>8}")
    print("-"*80)

    cat_results = []
    for cat, gene_list in func_categories.items():
        if not gene_list:
            continue
        # Concatenate all regions in this category
        combined = []
        for name, seq in sac_seqs:
            parts = []
            for g, s, e in gene_list:
                parts.append(seq[s:e+1])
            combined.append((name, ''.join(parts)))
        total_len = sum(e - s + 1 for _, s, e in gene_list)
        result = tajima_d(combined)
        if result:
            print(f"{cat:<20} {len(gene_list):>5} {total_len:>7} {result['S']:>6} "
                  f"{result['pi']:>10.4f} {result['theta_W']:>10.4f} {result['D']:>8.4f}")
            cat_results.append((cat, len(gene_list), total_len, result))

    # Also do the full locus
    result_full = tajima_d(sac_seqs)
    if result_full:
        print(f"\n{'FULL LOCUS':<20} {'all':>5} {len(sac_seqs[0][1]):>7} {result_full['S']:>6} "
              f"{result_full['pi']:>10.4f} {result_full['theta_W']:>10.4f} {result_full['D']:>8.4f}")

    # ── Analysis 3: Tajima's D per functional category, per species ──
    print("\n" + "="*80)
    print("TAJIMA'S D PER FUNCTIONAL CATEGORY — BY SPECIES")
    print("="*80)

    for sp in ['L_guanicoe', 'L_glama', 'V_vicugna', 'V_pacos']:
        sp_seqs = species_groups[sp]
        print(f"\n--- {sp} (n={len(sp_seqs)}) ---")
        print(f"{'Category':<20} {'Length':>7} {'S':>6} {'pi':>10} {'theta_W':>10} {'D':>8}")
        print("-"*70)
        for cat, gene_list in func_categories.items():
            if not gene_list:
                continue
            combined = []
            for name, seq in sp_seqs:
                parts = []
                for g, s, e in gene_list:
                    parts.append(seq[s:e+1])
                combined.append((name, ''.join(parts)))
            total_len = sum(e - s + 1 for _, s, e in gene_list)
            result = tajima_d(combined)
            if result and not math.isnan(result['D']):
                print(f"{cat:<20} {total_len:>7} {result['S']:>6} "
                      f"{result['pi']:>10.4f} {result['theta_W']:>10.4f} {result['D']:>8.4f}")
            elif result:
                print(f"{cat:<20} {total_len:>7} {result['S']:>6} "
                      f"{result['pi']:>10.4f} {result['theta_W']:>10.4f} {'NA':>8}")

        # Full locus per species
        result_sp = tajima_d(sp_seqs)
        if result_sp:
            D_str = f"{result_sp['D']:>8.4f}" if not math.isnan(result_sp['D']) else f"{'NA':>8}"
            print(f"{'FULL LOCUS':<20} {len(sp_seqs[0][1]):>7} {result_sp['S']:>6} "
                  f"{result_sp['pi']:>10.4f} {result_sp['theta_W']:>10.4f} {D_str}")

    # ── Analysis 4: Intergenic regions ──
    print("\n" + "="*80)
    print("TAJIMA'S D — INTERGENIC vs GENIC (ALL SAC)")
    print("="*80)

    # Build set of all genic positions
    genic_positions = set()
    for gene, start, end, gtype in regions:
        for pos in range(start, end+1):
            genic_positions.add(pos)

    total_len_aln = len(sac_seqs[0][1])
    intergenic_indices = [i for i in range(total_len_aln) if i not in genic_positions]
    genic_indices = sorted(genic_positions)

    # Extract intergenic
    intergenic_sub = []
    for name, seq in sac_seqs:
        intergenic_sub.append((name, ''.join(seq[i] for i in intergenic_indices)))

    genic_sub = []
    for name, seq in sac_seqs:
        genic_sub.append((name, ''.join(seq[i] for i in genic_indices)))

    print(f"\n{'Region':<20} {'Length':>7} {'S':>6} {'pi':>10} {'theta_W':>10} {'D':>8}")
    print("-"*70)

    r_inter = tajima_d(intergenic_sub)
    if r_inter:
        D_str = f"{r_inter['D']:>8.4f}" if not math.isnan(r_inter['D']) else "NA"
        print(f"{'Intergenic':<20} {len(intergenic_indices):>7} {r_inter['S']:>6} "
              f"{r_inter['pi']:>10.4f} {r_inter['theta_W']:>10.4f} {D_str:>8}")

    r_genic = tajima_d(genic_sub)
    if r_genic:
        D_str = f"{r_genic['D']:>8.4f}" if not math.isnan(r_genic['D']) else "NA"
        print(f"{'Genic (all)':<20} {len(genic_indices):>7} {r_genic['S']:>6} "
              f"{r_genic['pi']:>10.4f} {r_genic['theta_W']:>10.4f} {D_str:>8}")

    print("\nDone.")


if __name__ == '__main__':
    main()
