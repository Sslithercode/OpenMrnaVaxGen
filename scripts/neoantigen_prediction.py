
import os
# Force CPU-only: CUDA toolkit on this machine is missing libdevice.10.bc,
# which prevents TF JIT kernel compilation on the GPU. MHCflurry on CPU
# is fast enough for typical neoantigen candidate counts (<10k peptides).
os.environ["CUDA_VISIBLE_DEVICES"] = ""

import subprocess
import sys
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
from itertools import product

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).parent))
from paths import BASE, step_dir

VCF_FILE        = step_dir(2) / "filtered_variants.vcf.gz"
HLA_FILE        = step_dir(3) / "hla_alleles.txt"
REFERENCE       = BASE / "reference/hg38.fa"
OUT_DIR         = step_dir(4)
PEPTIDE_LENGTHS = [8, 9, 10, 11]
BINDING_CUTOFF  = 500  # nM IC50 — standard threshold for MHC binding

# ── Helpers ───────────────────────────────────────────────────────────────────

def load_hla_alleles(hla_file):
    """Load HLA alleles from Step 3 output."""
    with open(hla_file) as f:
        alleles = [line.strip() for line in f if line.strip()]
    print(f"Loaded {len(alleles)} HLA alleles: {alleles}")
    return alleles


def parse_vcf(vcf_file):
    """
    Parse VCF into a simple dataframe of mutations.
    Each row = one somatic variant with chrom, pos, ref, alt.
    """
    import pysam
    variants = []
    vcf = pysam.VariantFile(str(vcf_file))
    for rec in vcf.fetch():
        # Only keep PASS variants
        if rec.filter.keys() and "PASS" not in rec.filter.keys():
            continue
        # Only SNVs for now (ref and alt are single bases)
        ref = rec.ref
        alts = rec.alts
        if alts is None:
            continue
        for alt in alts:
            if len(ref) == 1 and len(alt) == 1:
                variants.append({
                    "chrom": rec.chrom,
                    "pos":   rec.pos,
                    "ref":   ref,
                    "alt":   alt,
                    "vaf":   _get_vaf(rec),
                })
    vcf.close()
    df = pd.DataFrame(variants)
    print(f"Loaded {len(df)} SNV variants from VCF")
    return df


def _get_vaf(rec):
    """Extract variant allele frequency from VCF record."""
    try:
        # Mutect2 stores AF in the tumor sample FORMAT field
        for sample in rec.samples.values():
            if "AF" in sample:
                return float(sample["AF"][0])
    except Exception:
        pass
    return None


def extract_peptides_from_vcf(variants_df, reference_fasta, peptide_lengths):
    """
    For each variant, extract the surrounding genomic context and
    generate all peptides of requested lengths that span the mutation.

    Uses samtools faidx to fetch reference sequence around each variant.
    Returns a dataframe of (mutation_id, peptide, is_mutant).
    """
    peptides = []
    max_len = max(peptide_lengths)

    for _, var in variants_df.iterrows():
        chrom = var["chrom"]
        pos   = var["pos"]  # 1-based in VCF
        ref   = var["ref"]
        alt   = var["alt"]

        # Fetch reference sequence window around the variant
        # We need max_len - 1 bases on each side to generate all peptides
        # DNA context needed: max_len * 3 codons worth
        flank = max_len * 3
        start = max(1, pos - flank)
        end   = pos + flank

        region = f"{chrom}:{start}-{end}"
        result = subprocess.run(
            ["samtools", "faidx", str(reference_fasta), region],
            capture_output=True, text=True
        )
        if result.returncode != 0:
            continue

        # Parse the FASTA sequence
        lines = result.stdout.strip().split("\n")
        if len(lines) < 2:
            continue
        ref_seq = "".join(lines[1:]).upper()

        # Find where our variant position falls in the window
        var_offset = pos - start  # 0-based offset in the window

        # Generate wildtype and mutant sequences
        wt_seq  = ref_seq
        mut_seq = ref_seq[:var_offset] + alt + ref_seq[var_offset + 1:]

        # Translate to amino acids and generate peptides
        for length in peptide_lengths:
            # Generate peptides that span the mutation position
            # In protein space: mutation falls at approximately var_offset // 3
            for frame in range(3):
                wt_aa  = _translate(wt_seq[frame:])
                mut_aa = _translate(mut_seq[frame:])

                if wt_aa == mut_aa:
                    continue  # synonymous mutation, skip

                # Find where the mutation falls in the AA sequence
                mut_pos_aa = (var_offset - frame) // 3
                if mut_pos_aa < 0:
                    continue

                # Generate all peptides of this length spanning the mutation
                for start_aa in range(
                    max(0, mut_pos_aa - length + 1),
                    min(mut_pos_aa + 1, len(mut_aa) - length + 1)
                ):
                    wt_pep  = wt_aa[start_aa:start_aa + length]
                    mut_pep = mut_aa[start_aa:start_aa + length]

                    if len(mut_pep) < length:
                        continue
                    if wt_pep == mut_pep:
                        continue
                    if "*" in mut_pep:
                        continue  # stop codon

                    mutation_id = f"{chrom}:{pos}:{ref}>{alt}"
                    peptides.append({
                        "mutation_id": mutation_id,
                        "peptide":     mut_pep,
                        "wt_peptide":  wt_pep,
                        "length":      length,
                        "vaf":         var["vaf"],
                    })

    df = pd.DataFrame(peptides).drop_duplicates(subset=["mutation_id", "peptide"])
    print(f"Generated {len(df)} candidate peptides")
    return df


def _translate(dna):
    """Translate DNA sequence to amino acids."""
    codon_table = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    }
    protein = ""
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i+3]
        aa = codon_table.get(codon, "X")
        if aa == "*":
            break
        protein += aa
    return protein


def predict_binding(peptides_df, hla_alleles):
    """
    Run MHCflurry presentation score predictions for all
    peptide/HLA combinations.
    """
    from mhcflurry import Class1PresentationPredictor
    predictor = Class1PresentationPredictor.load()

    print(f"Running MHCflurry on {len(peptides_df)} peptides x {len(hla_alleles)} alleles...")

    results = predictor.predict(
    peptides=peptides_df["peptide"].tolist(),
    alleles={"patient": hla_alleles},
    verbose=1,
    )

    # Merge predictions back with peptide metadata
    results = results.rename(columns={"peptide": "peptide"})
    merged = peptides_df.merge(results, on="peptide", how="left")

    print(f"Predictions complete. Shape: {merged.shape}")
    return merged


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load inputs
    hla_alleles = load_hla_alleles(HLA_FILE)
    variants_df = parse_vcf(VCF_FILE)

    if variants_df.empty:
        print("No variants found. Exiting.")
        sys.exit(1)

    # Generate candidate peptides
    peptides_df = extract_peptides_from_vcf(
        variants_df, REFERENCE, PEPTIDE_LENGTHS
    )

    if peptides_df.empty:
        print("No candidate peptides generated. Exiting.")
        sys.exit(1)

    # Predict MHC binding
    predictions_df = predict_binding(peptides_df, hla_alleles)

    # Save full results
    output_file = OUT_DIR / "candidate_neoantigens.tsv"
    predictions_df.to_csv(output_file, sep="\t", index=False)

    # Summary
    print(f"\n── Step 4 Summary ──")
    print(f"  Variants processed:      {len(variants_df)}")
    print(f"  Candidate peptides:      {len(peptides_df)}")
    print(f"  Predictions made:        {len(predictions_df)}")
    print(f"\nStep 4 complete. Output: {output_file}")


if __name__ == "__main__":
    main()