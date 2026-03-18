"""
Step 2: Somatic Variant Calling
Uses Mutect2 to find mutations present in tumor but not normal.
Output: a filtered VCF file containing candidate somatic variants.
"""

import subprocess
import sys
from pathlib import Path

GATK_JAR   = Path.home() / "melanoma-pipeline/tools/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar"
REFERENCE  = Path.home() / "melanoma-pipeline/reference/hg38.fa"
TUMOR_BAM  = Path.home() / "melanoma-pipeline/results/step1/sorted.bam"
NORMAL_BAM = Path.home() / "melanoma-pipeline/data/test/normal_chr17.bam"
OUT_DIR    = Path.home() / "melanoma-pipeline/results/step2"

def run(cmd, step_name):
    print(f"\n[{step_name}] running...")
    result = subprocess.run([str(c) for c in cmd])
    if result.returncode != 0:
        print(f"FAILED: {step_name}")
        sys.exit(1)
    print(f"DONE: {step_name}")

def gatk(args):
    return ["java", "-jar", str(GATK_JAR)] + args

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    raw_vcf      = OUT_DIR / "raw_variants.vcf.gz"
    filtered_vcf = OUT_DIR / "filtered_variants.vcf.gz"

    # Step 2a: Mutect2 — somatic variant caller
    # Compares tumor vs normal, flags positions where tumor differs
    run(gatk([
        "Mutect2",
        "-R",          str(REFERENCE),
        "-I",          str(TUMOR_BAM),
        "-I",          str(NORMAL_BAM),
        "-tumor",      "HCC1143_tumor",
        "-normal",     "HCC1143_normal",
        "-O",          str(raw_vcf),
        "--intervals", "chr17",
    ]), "Mutect2")

    # Step 2b: Filter low-confidence calls
    run(gatk([
        "FilterMutectCalls",
        "-R", str(REFERENCE),
        "-V", str(raw_vcf),
        "-O", str(filtered_vcf),
    ]), "FilterMutectCalls")

    # Count variants found
    result = subprocess.run(
        ["java", "-jar", str(GATK_JAR), "CountVariants", "-V", str(filtered_vcf)],
        capture_output=True, text=True
    )
    print(f"\n── Variant Count ──")
    print(result.stdout)
    print(f"Step 2 complete. Output: {filtered_vcf}")

if __name__ == "__main__":
    main()