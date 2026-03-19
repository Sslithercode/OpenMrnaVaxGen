"""
Step 1: Preprocessing
Marks duplicates and outputs a clean analysis-ready BAM.
"""

import subprocess
import sys
from pathlib import Path

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).parent))
from paths import BASE, step_dir, new_run

# ── Config ────────────────────────────────────────────────────────────────────

GATK_JAR  = BASE / "tools/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar"
REFERENCE = BASE / "reference/hg38.fa"
INPUT_BAM = BASE / "data/test/tumor_chr17.bam"
OUT_DIR   = step_dir(1)

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
    # Step 1 always starts a fresh run — downstream steps inherit the same run ID
    new_run()
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