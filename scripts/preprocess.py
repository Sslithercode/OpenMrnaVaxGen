"""
Step 1: Preprocessing
Marks duplicates and outputs a clean analysis-ready BAM.
"""

import subprocess
import sys
from pathlib import Path

# ── Config ─────────────────────────────────────────
# ───────────────────────────

GATK_JAR  = Path.home() / "melanoma-pipeline/tools/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar"
REFERENCE = Path.home() / "melanoma-pipeline/reference/b37.20.21.fasta"
INPUT_BAM = Path.home() / "melanoma-pipeline/data/test/test_tumor.bam"
OUT_DIR   = Path.home() / "melanoma-pipeline/results/step1"

# ── Helpers ───────────────────────────────────────────────────────────────────

def run(cmd, step_name):
    print(f"\n[{step_name}] running...")
    result = subprocess.run([str(c) for c in cmd])
    if result.returncode != 0:
        print(f"FAILED: {step_name}")
        sys.exit(1)
    print(f"DONE: {step_name}")

def gatk(args):
    return ["java", "-jar", str(GATK_JAR)] + args

# ── Pipeline ──────────────────────────────────────────────────────────────────

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Step 1a: Mark duplicates
    deduped_bam    = OUT_DIR / "marked_duplicates.bam"
    metrics_file   = OUT_DIR / "duplicate_metrics.txt"

    run(gatk([
        "MarkDuplicates",
        "--INPUT",        INPUT_BAM,
        "--OUTPUT",       deduped_bam,
        "--METRICS_FILE", metrics_file,
        "--VALIDATION_STRINGENCY", "SILENT",
        "--CREATE_INDEX", "true",
    ]), "MarkDuplicates")

    # Step 1b: Sort the BAM (required before BQSR)
    sorted_bam = OUT_DIR / "sorted.bam"

    run(gatk([
        "SortSam",
        "--INPUT",       deduped_bam,
        "--OUTPUT",      sorted_bam,
        "--SORT_ORDER",  "coordinate",
        "--CREATE_INDEX","true",
    ]), "SortSam")

    print(f"\nStep 1 complete. Output: {sorted_bam}")

if __name__ == "__main__":
    main()