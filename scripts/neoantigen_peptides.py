"""
Step 4a: Neoantigen Peptide Extraction
Parses somatic VCF and uses varcode + pyensembl to extract mutant peptides
from coding substitutions. Saves an intermediate TSV for Step 4b.

Kept in a separate process from Step 4b (MHCflurry) to avoid the
htslib / TensorFlow double-free crash on Linux.
"""

import os
import sys
import warnings
warnings.filterwarnings("ignore")

import pandas as pd

from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from paths import step_dir

VCF_FILE        = step_dir(2) / "filtered_variants.vcf.gz"
HLA_FILE        = step_dir(3) / "hla_alleles.txt"
OUT_DIR         = step_dir(4)
PEPTIDE_LENGTHS = [8, 9, 10, 11]


def load_hla_alleles(hla_file):
    with open(hla_file) as f:
        alleles = [line.strip() for line in f if line.strip()]
    print(f"Loaded {len(alleles)} HLA alleles: {alleles}")
    return alleles


def parse_vcf(vcf_file):
    """
    Pure Python gzip/text VCF parser — no pysam/htslib in this process.
    """
    import gzip
    variants = []
    opener = gzip.open if str(vcf_file).endswith(".gz") else open
    with opener(str(vcf_file), "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, pos, _, ref, alt_str, _, filter_, _, fmt = parts[:9]
            samples = parts[9:]
            if filter_ not in ("PASS", "."):
                continue
            format_keys = fmt.split(":")
            for i, alt in enumerate(alt_str.split(",")):
                if len(ref) != 1 or len(alt) != 1:
                    continue
                vaf = None
                for sample_data in samples:
                    d = dict(zip(format_keys, sample_data.split(":")))
                    if "AF" in d:
                        try:
                            af_values = d["AF"].split(",")
                            vaf = float(af_values[i]) if i < len(af_values) else None
                            break
                        except (ValueError, IndexError):
                            pass
                variants.append({
                    "chrom": chrom,
                    "pos":   int(pos),
                    "ref":   ref,
                    "alt":   alt,
                    "vaf":   vaf,
                })
    df = pd.DataFrame(variants)
    print(f"Loaded {len(df)} SNV variants from VCF")
    return df


def extract_peptides(variants_df, peptide_lengths):
    """
    Uses varcode + pyensembl (htslib-backed) to get true mutant protein
    sequences, then slices all peptides spanning the mutation position.
    """
    from pyensembl import EnsemblRelease
    from varcode import Variant
    from varcode.effects import Substitution

    genome = EnsemblRelease(109)
    peptides = []

    for _, var in variants_df.iterrows():
        chrom = str(var["chrom"])
        pos   = int(var["pos"])
        ref   = var["ref"]
        alt   = var["alt"]
        vaf   = var["vaf"]
        mutation_id = f"{chrom}:{pos}:{ref}>{alt}"

        ensembl_contig = chrom[3:] if chrom.startswith("chr") else chrom

        try:
            variant = Variant(
                contig=ensembl_contig, start=pos, ref=ref, alt=alt, genome=genome
            )
            effects = variant.effects()
            top = effects.top_priority_effect()
        except Exception:
            continue

        if not isinstance(top, Substitution):
            continue

        mut_protein = top.mutant_protein_sequence
        wt_protein  = top.original_protein_sequence
        if not mut_protein:
            continue

        mut_aa_pos = top.aa_mutation_start_offset
        gene_name  = top.gene.name

        for length in peptide_lengths:
            for start_aa in range(
                max(0, mut_aa_pos - length + 1),
                min(mut_aa_pos + 1, len(mut_protein) - length + 1)
            ):
                mut_pep = mut_protein[start_aa:start_aa + length]
                if len(mut_pep) < length or "*" in mut_pep:
                    continue
                wt_pep = wt_protein[start_aa:start_aa + length] if wt_protein else None
                if wt_pep == mut_pep:
                    continue
                peptides.append({
                    "mutation_id": mutation_id,
                    "peptide":     mut_pep,
                    "wt_peptide":  wt_pep,
                    "length":      length,
                    "vaf":         vaf,
                    "gene":        gene_name,
                })

    df = pd.DataFrame(peptides).drop_duplicates(subset=["mutation_id", "peptide", "wt_peptide"])
    print(f"Generated {len(df)} candidate peptides")
    return df


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    hla_alleles = load_hla_alleles(HLA_FILE)
    variants_df = parse_vcf(VCF_FILE)

    if variants_df.empty:
        print("No variants found. Exiting.")
        sys.exit(1)

    peptides_df = extract_peptides(variants_df, PEPTIDE_LENGTHS)

    if peptides_df.empty:
        print("No candidate peptides generated. Exiting.")
        sys.exit(1)

    # Save HLA alleles alongside peptides so Step 4b can read them
    hla_out = OUT_DIR / "hla_alleles_4b.txt"
    hla_out.write_text("\n".join(hla_alleles) + "\n")

    peptides_file = OUT_DIR / "candidate_peptides.tsv"
    peptides_df.to_csv(peptides_file, sep="\t", index=False)

    print(f"\n── Step 4a Summary ──")
    print(f"  Variants processed:  {len(variants_df)}")
    print(f"  Candidate peptides:  {len(peptides_df)}")
    print(f"\nStep 4a complete. Output: {peptides_file}")


if __name__ == "__main__":
    main()
