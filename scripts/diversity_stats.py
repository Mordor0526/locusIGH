#!/usr/bin/env python3
"""
Diversity statistics for the IGH locus — π ± SD, Hd ± SD, S, NH, Nsp, Tajima's D.
Formulae follow Nei (1987) and Tajima (1989).
Usage: python3 diversity_stats.py
Expects: IgH_103_aligned.msa in ~/locusigh/
"""

import os, math
from collections import defaultdict, Counter

ALN_FILE = os.path.expanduser("~/locusigh/IgH_103_aligned.msa")

def assign_species(name):
    n = name.lower()
    if 'cbac' in n or 'bactrianus' in n or n.startswith('cb'):
        return 'C_bactrianus'
    if n.endswith('g') and n[:-1].isdigit(): return 'L_guanicoe'
    if n.endswith('l') and n[:-1].isdigit(): return 'L_glama'
    if n.endswith('v') and n[:-1].isdigit(): return 'V_vicugna'
    if n.endswith('a') and n[:-1].isdigit(): return 'V_pacos'
    for tag, sp in [('guan','L_guanicoe'),('glam','L_glama'),('llam','L_glama'),
                    ('vicu','V_vicugna'),('paco','V_pacos'),('alpa','V_pacos')]:
        if tag in n: return sp
    return 'unknown'

def read_fasta(path):
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

def compute_diversity(sub_seqs):
    """
    Compute: n, L (usable sites), S, NH (num haplotypes), Nsp (num segregating patterns),
    pi ± SD, Hd ± SD, Tajima's D.
    """
    n = len(sub_seqs)
    if n < 2:
        return None
    seqlen = len(sub_seqs[0][1])

    # ── Per-site calculations ──
    S = 0
    usable_sites = 0
    total_pi = 0.0       # sum of per-site pairwise differences (normalised)
    total_pi_sq = 0.0    # for variance of pi

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
            counts = Counter(bases)
            site_diffs = 0
            allele_list = list(counts.keys())
            for a in range(len(allele_list)):
                for b in range(a+1, len(allele_list)):
                    site_diffs += counts[allele_list[a]] * counts[allele_list[b]]
            site_pairs = nn * (nn - 1) / 2
            pj = site_diffs / site_pairs
            total_pi += pj
            total_pi_sq += pj * pj

    # π per site
    pi = total_pi / usable_sites if usable_sites > 0 else 0

    # Variance of π — Nei (1987) eq. 10.7 (sampling variance)
    # Var(pi) = (1/L^2) * [ sum(pj^2) - (sum(pj))^2 / L ] * 2 / (n*(n-1))
    # Simplified: Nei & Li 1979, Nei 1987 ch.10
    # For nucleotide diversity, the SD reported by DnaSP uses:
    # SD(pi) = sqrt( (n+1)/(3*(n-1)*L) * pi  +  2*(n^2+n+3)/(9*n*(n-1)) * pi^2 / L )
    # This is based on Tajima (1983) and Nei (1987)
    if usable_sites > 0 and n > 1:
        L = usable_sites
        c1 = (n + 1) / (3 * (n - 1) * L)
        c2 = 2 * (n*n + n + 3) / (9 * n * (n - 1) * L)
        var_pi = c1 * pi + c2 * pi * pi
        sd_pi = math.sqrt(var_pi) if var_pi > 0 else 0
    else:
        sd_pi = 0

    # ── Haplotype diversity ──
    hap_counts = Counter()
    for name, seq in sub_seqs:
        # Use only usable sites for haplotype ID (but keep it simple: use full seq)
        hap_counts[seq] += 1
    NH = len(hap_counts)
    freqs = [c / n for c in hap_counts.values()]
    Hd = (n / (n - 1)) * (1 - sum(f*f for f in freqs))
    # SD of Hd — Nei (1987) eq. 8.4 (sampling variance)
    sum_fi2 = sum(f*f for f in freqs)
    sum_fi3 = sum(f*f*f for f in freqs)
    var_Hd = (2 / (n * (n-1))) * (
        2 * (n-2) * (sum_fi2**2 - sum(f**4 for f in freqs)) / ((n-2)*(n-3) if n > 3 else 1)
        - sum_fi2 * (1 - sum_fi2)
    ) if n > 3 else 0
    # Simpler Nei (1987) formula:
    # Var(Hd) = 2/(n(n-1)) * [ 2(n-2) * (Σfi³ - (Σfi²)²) + (Σfi²) - (Σfi²)² ]
    # But let's use the standard formula from Nei 1987 eq 8.12 (unbiased):
    var_Hd2 = (2 / (n*(n-1))) * (
        2*(n-2) * (sum(f**3 for f in freqs) - sum_fi2**2)
        + sum_fi2 - sum_fi2**2
    )
    sd_Hd = math.sqrt(abs(var_Hd2))

    # ── Tajima's D ──
    a1 = sum(1.0/i for i in range(1, n))
    a2 = sum(1.0/(i*i) for i in range(1, n))
    theta_W = S / a1 if a1 > 0 else 0

    D_val = float('nan')
    if S > 0 and n >= 4:
        b1 = (n + 1) / (3 * (n - 1))
        b2 = 2 * (n*n + n + 3) / (9 * n * (n - 1))
        c1_ = b1 - 1/a1
        c2_ = b2 - (n + 2) / (a1 * n) + a2 / (a1 * a1)
        e1 = c1_ / a1
        e2 = c2_ / (a1*a1 + a2)
        d = total_pi - theta_W
        var_d = e1 * S + e2 * S * (S - 1)
        if var_d > 0:
            D_val = d / math.sqrt(var_d)

    # Number of polymorphic sites (Nsp = S already)
    return {
        'n': n, 'L': usable_sites, 'S': S, 'NH': NH,
        'pi': pi, 'sd_pi': sd_pi,
        'Hd': Hd, 'sd_Hd': sd_Hd,
        'theta_W': theta_W / usable_sites if usable_sites > 0 else 0,
        'D': D_val
    }


def main():
    print("Reading alignment...")
    seqs = read_fasta(ALN_FILE)
    print(f"  {len(seqs)} sequences, {len(seqs[0][1])} sites\n")

    # Assign species
    groups = defaultdict(list)
    sac_seqs = []
    for name, seq in seqs:
        sp = assign_species(name)
        if sp == 'C_bactrianus':
            continue
        groups[sp].append((name, seq))
        sac_seqs.append((name, seq))

    # Ordered groups for output
    ordered = [
        ('All SAC', sac_seqs),
        ('L_guanicoe', groups['L_guanicoe']),
        ('L_glama', groups['L_glama']),
        ('V_vicugna', groups['V_vicugna']),
        ('V_pacos', groups['V_pacos']),
    ]

    print("="*100)
    print("DIVERSITY STATISTICS — FULL IGH LOCUS")
    print("="*100)
    print(f"\n{'Group':<16} {'n':>4} {'L':>8} {'S':>7} {'NH':>5} "
          f"{'pi':>10} {'SD(pi)':>10} {'Hd':>8} {'SD(Hd)':>8} "
          f"{'theta_W/site':>13} {'Tajima D':>9}")
    print("-"*100)

    for label, seqlist in ordered:
        r = compute_diversity(seqlist)
        if r:
            D_str = f"{r['D']:.4f}" if not math.isnan(r['D']) else 'NA'
            print(f"{label:<16} {r['n']:>4} {r['L']:>8} {r['S']:>7} {r['NH']:>5} "
                  f"{r['pi']:>10.6f} {r['sd_pi']:>10.6f} {r['Hd']:>8.4f} {r['sd_Hd']:>8.4f} "
                  f"{r['theta_W']:>13.6f} {D_str:>9}")

    # ── Also compute per-subspecies for VHH region ──
    # We need gene_coords_aln.tsv for VHH boundaries
    coords_file = os.path.expanduser("~/locusigh/gene_coords_aln.tsv")
    vhh_start, vhh_end = None, None
    if os.path.exists(coords_file):
        with open(coords_file) as f:
            header = f.readline()
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4 and 'vhh' in parts[0].lower():
                    vhh_start, vhh_end = int(parts[1]), int(parts[2])
                    break

    if vhh_start is not None:
        print("\n" + "="*100)
        print(f"DIVERSITY STATISTICS — VHH REGION (vhh3-1, positions {vhh_start}-{vhh_end})")
        print("="*100)
        print(f"\n{'Group':<20} {'n':>4} {'L':>6} {'S':>5} {'NH':>5} "
              f"{'pi':>10} {'SD(pi)':>10} {'Hd':>8} {'SD(Hd)':>8} {'Tajima D':>9}")
        print("-"*95)

        # Subspecies assignment for detailed VHH table
        def assign_subsp(name):
            n = name.lower()
            if 'cbac' in n or 'bactrianus' in n or n.startswith('cb'):
                return 'C_bactrianus'
            # Use the species-level first, then refine if possible
            sp = assign_species(name)
            # For subspecies, we'd need metadata — use species level
            return sp

        vhh_groups = defaultdict(list)
        vhh_all = []
        for name, seq in seqs:
            sp = assign_species(name)
            if sp == 'C_bactrianus':
                continue
            vhh_seq = seq[vhh_start:vhh_end+1]
            vhh_groups[sp].append((name, vhh_seq))
            vhh_all.append((name, vhh_seq))

        vhh_ordered = [
            ('All SAC', vhh_all),
            ('L_guanicoe', vhh_groups['L_guanicoe']),
            ('L_glama', vhh_groups['L_glama']),
            ('V_vicugna', vhh_groups['V_vicugna']),
            ('V_pacos', vhh_groups['V_pacos']),
        ]

        for label, seqlist in vhh_ordered:
            r = compute_diversity(seqlist)
            if r:
                D_str = f"{r['D']:.4f}" if not math.isnan(r['D']) else 'NA'
                print(f"{label:<20} {r['n']:>4} {r['L']:>6} {r['S']:>5} {r['NH']:>5} "
                      f"{r['pi']:>10.6f} {r['sd_pi']:>10.6f} {r['Hd']:>8.4f} {r['sd_Hd']:>8.4f} {D_str:>9}")

    print("\nDone.")


if __name__ == '__main__':
    main()
