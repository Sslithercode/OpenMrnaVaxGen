"""
Step 6: Epitope Ordering
Takes the top ranked neoantigens from Step 5 and finds the optimal
concatenation order to minimize junctional neoantigens.

When epitopes are joined together in a string, the junction between
peptide A and peptide B creates new subsequences that could:
- Bind MHC and trigger unwanted immune responses
- Create self-mimicking sequences

This is a shortest Hamiltonian path problem:
- Nodes = epitopes
- Edge weight = MHCflurry score of junctional peptide between A and B
- Goal = find ordering that minimizes total junction score

Between each epitope we insert a GPGPG linker which acts as a
spacer to further reduce junctional effects (BioNTech standard).

Output: ordered_epitopes.fasta — the epitope string ready for mRNA encoding
"""

import warnings
warnings.filterwarnings("ignore")

import sys
import pandas as pd
import numpy as np
import os

from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from paths import step_dir

# ── Config ────────────────────────────────────────────────────────────────────
os.environ["CUDA_VISIBLE_DEVICES"] = ""
STEP5_OUTPUT = step_dir(5) / "ranked_neoantigens.tsv"
HLA_FILE     = step_dir(3) / "hla_alleles.txt"
OUT_DIR      = step_dir(6)

LINKER       = "GPGPG"
JUNCTION_LEN = 9    # length of junctional peptides to check
MAX_EXACT    = 15   # use exact TSP for N <= 15, greedy above that.
                    # Held-Karp is O(2^N * N^2) — feasible up to ~15, unusable at 30.


# ── Junction scoring ──────────────────────────────────────────────────────────

def get_junction_peptides(pep_a, pep_b, length=JUNCTION_LEN):
    """
    Generate all peptides of given length that span the junction
    between peptide A, the GPGPG linker, and peptide B.

    Takes (length-1) residues from each side instead of hardcoded 4.
    Original used pep_a[-4:] + LINKER + pep_b[:4], which:
      - misses junctional peptides that draw more than 4 residues from one side
      - silently under-counts for 8-mer epitopes shorter than 4 residues on one side
    Correct window: (length-1) residues from tail of A + linker + (length-1) from head of B,
    which guarantees every possible length-mer spanning the junction is covered.
    """
    junction_seq = pep_a[-(length - 1):] + LINKER + pep_b[:(length - 1)]
    peptides = []
    for i in range(len(junction_seq) - length + 1):
        peptides.append(junction_seq[i:i + length])
    return peptides


def score_all_junctions(epitopes, hla_alleles):
    """
    For each ordered pair of epitopes (A, B), score the junctional peptides
    using MHCflurry. Returns an N x N matrix of junction scores.
    Lower = better (less likely to create unwanted neoantigens).
    junction_matrix[i][j] = worst junctional score if epitope i is followed by epitope j.
    """
    from mhcflurry import Class1PresentationPredictor
    predictor = Class1PresentationPredictor.load()

    n = len(epitopes)
    junction_matrix = np.zeros((n, n))

    print(f"Scoring junctions for {n} epitopes ({n * (n-1)} ordered pairs)...")

    # Collect all unique junction peptides across all pairs in one pass
    junction_map = {}
    all_junction_peptides = set()

    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            junctions = get_junction_peptides(epitopes[i], epitopes[j])
            junction_map[(i, j)] = junctions
            all_junction_peptides.update(junctions)

    all_junction_peptides = list(all_junction_peptides)

    if not all_junction_peptides:
        print("No junction peptides generated.")
        return junction_matrix

    # Score all junction peptides in a single batch call — avoids N^2 model loads
    results = predictor.predict(
        peptides=all_junction_peptides,
        alleles={"patient": hla_alleles},
        verbose=0,
    )

    # Build lookup: peptide -> worst (max) presentation score across alleles
    pep_scores = results.groupby("peptide")["presentation_score"].max().to_dict()

    # Fill matrix: each cell = worst junctional peptide score for that pair
    for (i, j), junctions in junction_map.items():
        scores = [pep_scores.get(p, 0) for p in junctions]
        junction_matrix[i][j] = max(scores) if scores else 0

    print(f"Junction scoring complete.")
    print(f"  Max junction score:  {junction_matrix.max():.3f}")
    print(f"  Mean junction score: {junction_matrix[junction_matrix > 0].mean():.3f}")
    return junction_matrix


# ── Ordering algorithms ───────────────────────────────────────────────────────

def greedy_ordering(epitopes, junction_matrix):
    """
    Greedy nearest-neighbor TSP approximation.
    At each step, pick the next epitope that creates the lowest junction
    score with the current last epitope.
    Used when N > MAX_EXACT (typically N > 15).
    """
    n = len(epitopes)
    unvisited = set(range(n))

    # Start from the epitope with the lowest average outgoing junction score
    avg_out = [junction_matrix[i].sum() for i in range(n)]
    start = int(np.argmin(avg_out))

    path = [start]
    unvisited.remove(start)

    while unvisited:
        current = path[-1]
        best_next = min(unvisited, key=lambda j: junction_matrix[current][j])
        path.append(best_next)
        unvisited.remove(best_next)

    total_score = sum(
        junction_matrix[path[i]][path[i + 1]]
        for i in range(len(path) - 1)
    )
    return path, total_score


def exact_ordering(epitopes, junction_matrix):
    """
    Exact TSP solution via Held-Karp dynamic programming.
    O(2^N * N^2) time — only feasible for N <= MAX_EXACT (~15).
    Finds the globally optimal epitope ordering.

    dp[mask][i] = minimum total junction score to visit exactly the epitopes
                  in `mask`, ending at epitope i.
    """
    n = len(epitopes)

    print(f"Running exact Held-Karp TSP for N={n} epitopes...")

    INF  = float('inf')
    dp   = [[INF] * n for _ in range(1 << n)]
    prev = [[-1]  * n for _ in range(1 << n)]

    # Base case: start at any single epitope with zero cost
    for i in range(n):
        dp[1 << i][i] = 0

    for mask in range(1, 1 << n):
        for last in range(n):
            if not (mask & (1 << last)):
                continue
            if dp[mask][last] == INF:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                new_mask = mask | (1 << nxt)
                new_cost = dp[mask][last] + junction_matrix[last][nxt]
                if new_cost < dp[new_mask][nxt]:
                    dp[new_mask][nxt]   = new_cost
                    prev[new_mask][nxt] = last

    # Find best terminal epitope across all complete paths
    full_mask  = (1 << n) - 1
    best_last  = min(range(n), key=lambda i: dp[full_mask][i])
    best_score = dp[full_mask][best_last]

    # Reconstruct path by backtracking through prev pointers
    path    = []
    mask    = full_mask
    current = best_last
    while current != -1:
        path.append(current)
        last    = prev[mask][current]
        mask   ^= (1 << current)
        current = last
    path.reverse()

    return path, best_score


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(STEP5_OUTPUT, sep="\t")
    print(f"Loaded {len(df)} ranked epitopes from Step 5")

    with open(HLA_FILE) as f:
        hla_alleles = [l.strip() for l in f if l.strip()]
    print(f"HLA alleles: {hla_alleles}")

    # Deduplicate by mutation_id — keep best-ranked peptide per mutation.
    # Step 5 can surface multiple peptide lengths from the same somatic mutation.
    # Including both doubles antigen load with no additional epitope diversity.
    n_before = len(df)
    df = df.groupby("mutation_id", sort=False).first().reset_index()
    df = df.sort_values("composite_score", ascending=False).reset_index(drop=True)
    n_after = len(df)
    if n_before != n_after:
        print(f"Deduplicated by mutation_id: {n_before} -> {n_after} epitopes "
              f"({n_before - n_after} redundant lengths removed)")

    epitopes = df["peptide"].tolist()
    n = len(epitopes)

    if n == 0:
        print("No epitopes to order. Exiting.")
        return

    # Score all junctions
    junction_matrix = score_all_junctions(epitopes, hla_alleles)

    # Choose algorithm: exact for small N, greedy for large N
    print(f"\nFinding optimal epitope ordering (N={n})...")
    if n <= MAX_EXACT:
        path, total_score = exact_ordering(epitopes, junction_matrix)
        method = "exact TSP (Held-Karp)"
    else:
        path, total_score = greedy_ordering(epitopes, junction_matrix)
        method = "greedy nearest-neighbor"

    print(f"Ordering method: {method}")
    print(f"Total junction score: {total_score:.4f} (lower = better)")

    # Build ordered sequence
    ordered_epitopes = [epitopes[i] for i in path]
    ordered_df = df.iloc[path].copy().reset_index(drop=True)
    ordered_df["order"] = range(1, len(path) + 1)

    # Compute per-junction scores along the chosen path
    junction_scores_path = [
        junction_matrix[path[i]][path[i + 1]]
        for i in range(len(path) - 1)
    ]

    print(f"\n── Ordered Epitopes ──")
    for idx, ep in enumerate(ordered_epitopes):
        allele = ordered_df.iloc[idx]["best_allele"]
        flag   = str(ordered_df.iloc[idx].get("flags", ""))
        flag_s = f"  [{flag}]" if flag and flag != "nan" else ""
        if idx < len(junction_scores_path):
            jsc = junction_scores_path[idx]
            print(f"  {idx+1:2d}. {ep:12s}  ({allele})  -> junction: {jsc:.3f}{flag_s}")
        else:
            print(f"  {idx+1:2d}. {ep:12s}  ({allele}){flag_s}")

    # Join with GPGPG linkers
    epitope_string = LINKER.join(ordered_epitopes)
    print(f"\n── Epitope String ──")
    print(f"  Epitopes: {len(ordered_epitopes)}")
    print(f"  Linker:   {LINKER}")
    print(f"  Length:   {len(epitope_string)} aa")
    print(f"\n  {epitope_string}")

    # Save outputs
    ordered_df.to_csv(OUT_DIR / "ordered_epitopes.tsv", sep="\t", index=False)

    fasta_file = OUT_DIR / "ordered_epitopes.fasta"
    with open(fasta_file, "w") as f:
        f.write(
            f">vaccine_epitope_string | {len(ordered_epitopes)} epitopes | "
            f"linker={LINKER} | method={method} | junction_score={total_score:.4f}\n"
        )
        for i in range(0, len(epitope_string), 60):
            f.write(epitope_string[i:i + 60] + "\n")

    # Save junction matrix for inspection / Step 9 report
    junction_df = pd.DataFrame(
        junction_matrix,
        index=epitopes,
        columns=epitopes,
    )
    junction_df.to_csv(OUT_DIR / "junction_scores.tsv", sep="\t")

    print(f"\n── Step 6 Summary ──")
    print(f"  Input epitopes:    {n_before}")
    print(f"  After dedup:       {n_after}")
    print(f"  Ordering method:   {method}")
    print(f"  Total junc. score: {total_score:.4f}")
    if junction_scores_path:
        worst_pos = junction_scores_path.index(max(junction_scores_path))
        print(f"  Worst junction:    {max(junction_scores_path):.3f}  "
              f"(position {worst_pos + 1} -> {worst_pos + 2}: "
              f"{ordered_epitopes[worst_pos]} | {ordered_epitopes[worst_pos + 1]})")
    else:
        print(f"  Worst junction:    N/A (single epitope)")
    print(f"  Epitope string:    {len(epitope_string)} aa")
    print(f"\nStep 6 complete.")
    print(f"  Ordered TSV: {OUT_DIR / 'ordered_epitopes.tsv'}")
    print(f"  FASTA:       {fasta_file}")
    print(f"  Junctions:   {OUT_DIR / 'junction_scores.tsv'}")


if __name__ == "__main__":
    main()