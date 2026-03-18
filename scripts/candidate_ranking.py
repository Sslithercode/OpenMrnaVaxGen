"""
Step 5: Immunogenicity Scoring & Ranking (Literature-informed)
Takes MHCflurry predictions from Step 4 and ranks candidates using
a composite score based on 2024 literature consensus (IMPROVE, ImmuneMirror):

  1. Presentation score (MHCflurry)     — most important per IMPROVE paper
  2. Agretopicity (DAI)                 — mutant/wildtype affinity ratio, high importance
  3. VAF/clonality                      — high importance per CheckMate 153 data
  4. BLOSUM mutation score              — medium importance
  5. Foreignness (BLOSUM self-similarity)— minor but adds signal

References:
- IMPROVE (Frontiers Immunology 2024): feature importance analysis
- ImmuneMirror (Briefings Bioinformatics 2024): AUC 0.64 → 0.87 with these features
- NeoPrecis (bioRxiv 2025): agretopicity + clonality are strongest TCR recognition proxies
- BioNTech/Moderna clinical standard: VAF >= 0.05 hard filter (clonal mutation requirement)
"""

import warnings, os
warnings.filterwarnings("ignore")
os.environ["CUDA_VISIBLE_DEVICES"] = ""
import pandas as pd
import numpy as np

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).parent))
from paths import step_dir

# ── Config ────────────────────────────────────────────────────────────────────

STEP4_OUTPUT = step_dir(4) / "candidate_neoantigens.tsv"
HLA_FILE     = step_dir(3) / "hla_alleles.txt"
OUT_DIR      = step_dir(5)

# Weights based on IMPROVE feature importance analysis
W_PRESENTATION = 0.40
W_AGRETOPICITY = 0.25
W_VAF          = 0.20
W_BLOSUM       = 0.10
W_FOREIGNNESS  = 0.05

# Filters
MAX_AFFINITY_NM  = 500
MIN_PRESENTATION = 0.1
MIN_VAF          = 0.05   # [ADDED] Clinical standard: BioNTech/Moderna filter subclonal mutations.
                           # VAF < 0.05 = mutation present in <5% of tumor cells -> poor vaccine target.
TOP_N            = 30

# Flag thresholds (candidates kept but annotated for wet lab review)
FLAG_NEG_AGRETOPICITY = 0.0   # [ADDED] Negative agretopicity = mutant binds MHC *worse* than wildtype.
                               # T cell repertoire is not depleted of wildtype-specific cells,
                               # so immune response may be misdirected. Flag for review.

# BLOSUM62 matrix
BLOSUM62 = {
    ('A','A'):4,  ('A','R'):-1, ('A','N'):-2, ('A','D'):-2, ('A','C'):0,
    ('A','Q'):-1, ('A','E'):-1, ('A','G'):0,  ('A','H'):-2, ('A','I'):-1,
    ('A','L'):-1, ('A','K'):-1, ('A','M'):-1, ('A','F'):-2, ('A','P'):-1,
    ('A','S'):1,  ('A','T'):0,  ('A','W'):-3, ('A','Y'):-2, ('A','V'):0,
    ('R','R'):5,  ('R','N'):-1, ('R','D'):-2, ('R','C'):-3, ('R','Q'):1,
    ('R','E'):0,  ('R','G'):-2, ('R','H'):0,  ('R','I'):-3, ('R','L'):-2,
    ('R','K'):2,  ('R','M'):-1, ('R','F'):-3, ('R','P'):-2, ('R','S'):-1,
    ('R','T'):-1, ('R','W'):-3, ('R','Y'):-2, ('R','V'):-3,
    ('N','N'):6,  ('N','D'):1,  ('N','C'):-3, ('N','Q'):0,  ('N','E'):0,
    ('N','G'):0,  ('N','H'):1,  ('N','I'):-3, ('N','L'):-3, ('N','K'):0,
    ('N','M'):-2, ('N','F'):-3, ('N','P'):-2, ('N','S'):1,  ('N','T'):0,
    ('N','W'):-4, ('N','Y'):-2, ('N','V'):-3,
    ('D','D'):6,  ('D','C'):-3, ('D','Q'):0,  ('D','E'):2,  ('D','G'):-1,
    ('D','H'):-1, ('D','I'):-3, ('D','L'):-4, ('D','K'):-1, ('D','M'):-3,
    ('D','F'):-3, ('D','P'):-1, ('D','S'):0,  ('D','T'):-1, ('D','W'):-4,
    ('D','Y'):-3, ('D','V'):-3,
    ('C','C'):9,  ('C','Q'):-3, ('C','E'):-4, ('C','G'):-3, ('C','H'):-3,
    ('C','I'):-1, ('C','L'):-1, ('C','K'):-3, ('C','M'):-1, ('C','F'):-2,
    ('C','P'):-3, ('C','S'):-1, ('C','T'):-1, ('C','W'):-2, ('C','Y'):-2,
    ('C','V'):-1,
    ('Q','Q'):5,  ('Q','E'):2,  ('Q','G'):-2, ('Q','H'):0,  ('Q','I'):-3,
    ('Q','L'):-2, ('Q','K'):1,  ('Q','M'):0,  ('Q','F'):-3, ('Q','P'):-1,
    ('Q','S'):0,  ('Q','T'):-1, ('Q','W'):-2, ('Q','Y'):-1, ('Q','V'):-2,
    ('E','E'):5,  ('E','G'):-2, ('E','H'):0,  ('E','I'):-3, ('E','L'):-3,
    ('E','K'):1,  ('E','M'):-2, ('E','F'):-3, ('E','P'):-1, ('E','S'):0,
    ('E','T'):-1, ('E','W'):-3, ('E','Y'):-2, ('E','V'):-2,
    ('G','G'):6,  ('G','H'):-2, ('G','I'):-4, ('G','L'):-4, ('G','K'):-2,
    ('G','M'):-3, ('G','F'):-3, ('G','P'):-2, ('G','S'):0,  ('G','T'):-2,
    ('G','W'):-2, ('G','Y'):-3, ('G','V'):-3,
    ('H','H'):8,  ('H','I'):-3, ('H','L'):-3, ('H','K'):-1, ('H','M'):-2,
    ('H','F'):-1, ('H','P'):-2, ('H','S'):-1, ('H','T'):-2, ('H','W'):-2,
    ('H','Y'):2,  ('H','V'):-3,
    ('I','I'):4,  ('I','L'):2,  ('I','K'):-1, ('I','M'):1,  ('I','F'):0,
    ('I','P'):-3, ('I','S'):-2, ('I','T'):-1, ('I','W'):-3, ('I','Y'):-1,
    ('I','V'):3,
    ('L','L'):4,  ('L','K'):-2, ('L','M'):2,  ('L','F'):0,  ('L','P'):-3,
    ('L','S'):-2, ('L','T'):-1, ('L','W'):-2, ('L','Y'):-1, ('L','V'):1,
    ('K','K'):5,  ('K','M'):-1, ('K','F'):-3, ('K','P'):-1, ('K','S'):0,
    ('K','T'):-1, ('K','W'):-3, ('K','Y'):-2, ('K','V'):-2,
    ('M','M'):5,  ('M','F'):0,  ('M','P'):-2, ('M','S'):-1, ('M','T'):-1,
    ('M','W'):-1, ('M','Y'):-1, ('M','V'):1,
    ('F','F'):6,  ('F','P'):-4, ('F','S'):-2, ('F','T'):-2, ('F','W'):1,
    ('F','Y'):3,  ('F','V'):-1,
    ('P','P'):7,  ('P','S'):-1, ('P','T'):-1, ('P','W'):-4, ('P','Y'):-3,
    ('P','V'):-2,
    ('S','S'):4,  ('S','T'):1,  ('S','W'):-3, ('S','Y'):-2, ('S','V'):-2,
    ('T','T'):5,  ('T','W'):-2, ('T','Y'):-2, ('T','V'):0,
    ('W','W'):11, ('W','Y'):2,  ('W','V'):-3,
    ('Y','Y'):7,  ('Y','V'):-1,
    ('V','V'):4,
}

def blosum62(aa1, aa2):
    key = (aa1, aa2)
    if key in BLOSUM62:
        return BLOSUM62[key]
    return BLOSUM62.get((aa2, aa1), 0)


# ── Feature Computation ───────────────────────────────────────────────────────

def compute_agretopicity(df, hla_alleles):
    """
    Agretopicity (DAI) = log2(wt_affinity / mut_affinity)
    Higher = mutant binds MHC much better than wildtype = better vaccine target.
    Wildtype-specific T cells were deleted during thymic selection,
    so only T cells recognizing the mutant survive in the repertoire.
    Negative values: mutant binds MHC *worse* than wildtype — flagged for review.
    """
    print("Computing agretopicity (running MHCflurry on wildtype peptides)...")
    from mhcflurry import Class1PresentationPredictor
    predictor = Class1PresentationPredictor.load()

    valid_wt = df[
        df["wt_peptide"].str.len() == df["length"]
    ]["wt_peptide"].dropna().unique().tolist()

    if not valid_wt:
        df["agretopicity"] = 0.0
        return df

    wt_results = predictor.predict(
        peptides=valid_wt,
        alleles={"patient": hla_alleles},
        verbose=0,
    )
    wt_best = wt_results.groupby("peptide")["affinity"].min().reset_index()
    wt_best.columns = ["wt_peptide", "wt_affinity"]

    df = df.merge(wt_best, on="wt_peptide", how="left")
    df["wt_affinity"] = df["wt_affinity"].fillna(50000)
    df["agretopicity"] = np.log2(
        df["wt_affinity"].clip(1) / df["affinity"].clip(1)
    )
    print(f"  Mean agretopicity: {df['agretopicity'].mean():.2f}")
    return df


def compute_blosum_score(df):
    """BLOSUM62 score at mutated positions. Lower = more radical substitution."""
    scores = []
    for _, row in df.iterrows():
        mut, wt = row["peptide"], row["wt_peptide"]
        if len(mut) != len(wt):
            scores.append(0)
            continue
        changed = [blosum62(a, b) for a, b in zip(wt, mut) if a != b]
        scores.append(np.mean(changed) if changed else 0)
    df["blosum_score"] = scores
    return df


def compute_foreignness(df):
    """1 - BLOSUM kernel similarity between mutant and wildtype peptide."""
    def kernel_sim(p1, p2):
        if len(p1) != len(p2):
            return 0.0
        score     = sum(blosum62(a, b) for a, b in zip(p1, p2))
        self_score = sum(blosum62(a, a) for a in p1)
        return score / (self_score + 1e-9)

    df["foreignness"] = df.apply(
        lambda r: max(0, 1 - kernel_sim(r["peptide"], r["wt_peptide"])), axis=1
    )
    return df


# [ADDED] Flag candidates for wet lab review without removing them from the output.
# Flags are informational — final include/exclude decision belongs to the researcher.
def add_flags(df):
    """
    Annotate candidates with QC flags. Flagged candidates are kept in the output
    but marked for manual review before synthesis.

    Flags:
      LOW_VAF          — VAF between MIN_VAF and 0.10. Passed the hard filter but
                         still moderately subclonal. T cell coverage may be incomplete.
      NEG_AGRETOPICITY — Agretopicity < 0. Mutant binds MHC worse than wildtype.
                         Wildtype-specific T cells not depleted -> tolerance risk.
    """
    flags = []
    for _, row in df.iterrows():
        f = []
        vaf = row.get("vaf", None)
        if vaf is not None and not pd.isna(vaf) and vaf < 0.10:
            f.append("LOW_VAF")
        if row.get("agretopicity", 0) < FLAG_NEG_AGRETOPICITY:
            f.append("NEG_AGRETOPICITY")
        flags.append("|".join(f) if f else "")
    df["flags"] = flags
    return df


def score_candidates(df):
    df["score_presentation"] = df["presentation_score"].clip(0, 1)

    ag = df["agretopicity"]
    df["score_agretopicity"] = (ag - ag.min()) / (ag.max() - ag.min() + 1e-9)

    vaf = df["vaf"].fillna(0.5)
    df["score_vaf"] = (vaf - vaf.min()) / (vaf.max() - vaf.min() + 1e-9)

    bl = df["blosum_score"]
    df["score_blosum"] = 1 - (bl - bl.min()) / (bl.max() - bl.min() + 1e-9)

    df["score_foreignness"] = df["foreignness"].clip(0, 1)

    df["composite_score"] = (
        W_PRESENTATION * df["score_presentation"] +
        W_AGRETOPICITY * df["score_agretopicity"] +
        W_VAF          * df["score_vaf"] +
        W_BLOSUM       * df["score_blosum"] +
        W_FOREIGNNESS  * df["score_foreignness"]
    )
    return df


def filter_candidates(df):
    n = len(df)
    df = df[df["affinity"] <= MAX_AFFINITY_NM]
    print(f"  After affinity filter (<{MAX_AFFINITY_NM} nM): {len(df)} / {n}")
    df = df[df["presentation_score"] >= MIN_PRESENTATION]
    print(f"  After presentation filter (>{MIN_PRESENTATION}): {len(df)}")

    # [ADDED] Hard VAF filter — clinical standard (BioNTech/Moderna pipelines).
    # Removes subclonal mutations where <5% of tumor cells carry the variant.
    # These would train T cells against targets absent in the vast majority of tumor.
    # Candidates failing this filter are saved separately for transparency/audit trail.
    n_before_vaf = len(df)
    vaf_present = df["vaf"].notna()
    low_vaf_mask = vaf_present & (df["vaf"] < MIN_VAF)
    df_low_vaf = df[low_vaf_mask].copy()
    df = df[~low_vaf_mask]
    print(f"  After VAF filter (>={MIN_VAF}): {len(df)} / {n_before_vaf} "
          f"({len(df_low_vaf)} removed as subclonal)")

    # [CHANGED] Returns tuple (df, df_low_vaf) — caller unpacks both
    return df, df_low_vaf


def select_top_n(df, top_n):
    df = df.sort_values("composite_score", ascending=False)
    best = df.groupby("mutation_id").first().reset_index()
    best = best.sort_values("composite_score", ascending=False)
    if len(best) >= top_n:
        return best.head(top_n)
    filler = df[~df.index.isin(best.index)].head(top_n - len(best))
    return pd.concat([best, filler]).sort_values(
        "composite_score", ascending=False
    ).head(top_n)


# [ADDED] HLA coverage report — how many top candidates hit each allele.
# A balanced vaccine should cover all patient alleles, not just the most common.
# If one allele has 0 candidates in the top N, the patient's T cells expressing
# that allele have no vaccine targets — may need to relax filters for that allele.
def print_hla_coverage(df, hla_alleles):
    print("\n── HLA Coverage (top candidates) ──")
    for allele in hla_alleles:
        count = (df["best_allele"] == allele).sum()
        bar = "█" * count
        print(f"  {allele:20s}: {count:3d}  {bar}")
    alleles_with_zero = [a for a in hla_alleles if (df["best_allele"] == a).sum() == 0]
    if alleles_with_zero:
        print(f"\n  WARNING: No candidates for: {', '.join(alleles_with_zero)}")
        print(f"     Consider relaxing affinity/presentation thresholds for these alleles.")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(STEP4_OUTPUT, sep="\t")
    print(f"Loaded {len(df)} candidates from Step 4")

    with open(HLA_FILE) as f:
        hla_alleles = [l.strip() for l in f if l.strip()]
    print(f"HLA alleles: {hla_alleles}")

    print("\nApplying filters:")
    # [CHANGED] filter_candidates now returns (df, df_low_vaf) — two values
    df, df_low_vaf = filter_candidates(df)
    if df.empty:
        print("No candidates passed filters.")
        return

    print("\nComputing features:")
    df = compute_agretopicity(df, hla_alleles)
    df = compute_blosum_score(df)
    df = compute_foreignness(df)
    df = score_candidates(df)

    # [ADDED] Flag candidates before selecting top N
    df = add_flags(df)

    top = select_top_n(df, TOP_N)

    # Save outputs
    df.sort_values("composite_score", ascending=False).to_csv(
        OUT_DIR / "all_scored_candidates.tsv", sep="\t", index=False
    )
    top.to_csv(OUT_DIR / "ranked_neoantigens.tsv", sep="\t", index=False)

    # [ADDED] Save subclonal candidates separately for transparency/audit trail
    if not df_low_vaf.empty:
        df_low_vaf.to_csv(
            OUT_DIR / "subclonal_filtered.tsv", sep="\t", index=False
        )
        print(f"\n  Subclonal candidates (VAF < {MIN_VAF}) saved to: "
              f"{OUT_DIR / 'subclonal_filtered.tsv'}")

    # [CHANGED] Added flags column to display
    print(f"\n── Top {TOP_N} Vaccine Candidates ──")
    cols = ["mutation_id", "peptide", "best_allele",
            "affinity", "presentation_score", "agretopicity",
            "vaf", "composite_score", "flags"]
    print(top[cols].to_string(index=False))

    # [ADDED] Flag summary — surfaces issues without burying them in the table
    flagged = top[top["flags"] != ""]
    if not flagged.empty:
        print(f"\n── Flagged Candidates ({len(flagged)}) ──")
        print("  These passed all filters but warrant manual review before synthesis:")
        for _, row in flagged.iterrows():
            print(f"  {row['peptide']:12s}  {row['flags']}")

    # [ADDED] HLA coverage report
    print_hla_coverage(top, hla_alleles)

    print(f"\n── Step 5 Summary ──")
    print(f"  Input:              {len(pd.read_csv(STEP4_OUTPUT, sep=chr(9)))}")
    print(f"  After filters:      {len(df)}")
    print(f"  Subclonal removed:  {len(df_low_vaf)} (VAF < {MIN_VAF})")
    flagged_count = (top["flags"] != "").sum()
    print(f"  Flagged in top {TOP_N}:  {flagged_count}")
    print(f"  Vaccine targets:    {len(top)}")
    print(f"  Unique mutations:   {top['mutation_id'].nunique()}")
    print(f"\nStep 5 complete.")
    print(f"  Ranked:     {OUT_DIR / 'ranked_neoantigens.tsv'}")
    print(f"  All scored: {OUT_DIR / 'all_scored_candidates.tsv'}")
    if not df_low_vaf.empty:
        print(f"  Subclonal:  {OUT_DIR / 'subclonal_filtered.tsv'}")


if __name__ == "__main__":
    main()