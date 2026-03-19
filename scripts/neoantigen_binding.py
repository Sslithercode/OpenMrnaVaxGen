"""
Step 4b: MHC Binding Prediction
Loads candidate peptides from Step 4a and runs MHCflurry presentation
score predictions. Kept in a separate process from Step 4a (varcode/pyensembl)
to avoid the htslib / TensorFlow double-free crash on Linux.
"""

import os
os.environ["CUDA_VISIBLE_DEVICES"] = ""

import sys
import warnings
warnings.filterwarnings("ignore")

import pandas as pd

from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from paths import step_dir

OUT_DIR        = step_dir(4)
PEPTIDES_FILE  = OUT_DIR / "candidate_peptides.tsv"
HLA_FILE       = OUT_DIR / "hla_alleles_4b.txt"
BINDING_CUTOFF = 500  # nM IC50 — applied here; Step 5 also filters at this threshold


def main():
    if not PEPTIDES_FILE.exists():
        print(f"Peptides file not found: {PEPTIDES_FILE}")
        print("Run Step 4a first.")
        sys.exit(1)

    peptides_df = pd.read_csv(PEPTIDES_FILE, sep="\t")
    print(f"Loaded {len(peptides_df)} candidate peptides")

    with open(HLA_FILE) as f:
        hla_alleles = [l.strip() for l in f if l.strip()]
    print(f"HLA alleles: {hla_alleles}")

    from mhcflurry import Class1PresentationPredictor
    predictor = Class1PresentationPredictor.load()

    print(f"Running MHCflurry on {len(peptides_df)} peptides x {len(hla_alleles)} alleles...")
    results = predictor.predict(
        peptides=peptides_df["peptide"].tolist(),
        alleles={"patient": hla_alleles},
        verbose=1,
    )

    merged = peptides_df.merge(results, on="peptide", how="left")

    before = len(merged)
    merged = merged[merged["affinity"] <= BINDING_CUTOFF].reset_index(drop=True)
    print(f"After {BINDING_CUTOFF} nM cutoff: {len(merged)} / {before} peptides retained")

    output_file = OUT_DIR / "candidate_neoantigens.tsv"
    merged.to_csv(output_file, sep="\t", index=False)

    print(f"\n── Step 4b Summary ──")
    print(f"  Peptides input:      {before}")
    print(f"  Binders retained:    {len(merged)}")
    print(f"\nStep 4b complete. Output: {output_file}")


if __name__ == "__main__":
    main()