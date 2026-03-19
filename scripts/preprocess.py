"""
Step 1: Preprocessing
Sorts input BAM, marks duplicates, and outputs a clean analysis-ready BAM.
"""
 
import subprocess
import sys
from pathlib import Path
 
sys.path.insert(0, str(Path(__file__).parent))
from paths import BASE, step_dir, new_run
 
# ── Config ────────────────────────────────────────────────────────────────────
 
GATK_JAR  = BASE / "tools/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar"
REFERENCE = BASE / "reference/hg38.fa"
INPUT_BAM = BASE / "data/test/tumor_chr17.bam"
OUT_DIR   = step_dir(1)
 
# ── Helpers ───────────────────────────────────────────────────────────────────
 
def run(cmd, step_name):
    print(f"\n[{step_name}] running...")
    result = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    if result.returncode != 0:
        print(f"FAILED: {step_name}")
        print(result.stderr)
        sys.exit(1)
    print(f"DONE: {step_name}")
 
def gatk(args):
    return ["java", "-jar", str(GATK_JAR)] + args
 
# ── Pipeline ──────────────────────────────────────────────────────────────────
 
def main():
    # Step 1 always starts a fresh run — downstream steps inherit the same run ID
    new_run()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
 
    # Step 1a: Sort by coordinate first — required for correct duplicate marking
    sorted_bam = OUT_DIR / "sorted.bam"
 
    run(gatk([
        "SortSam",
        "--INPUT",        INPUT_BAM,
        "--OUTPUT",       sorted_bam,
        "--SORT_ORDER",   "coordinate",
        "--CREATE_INDEX", "true",
    ]), "SortSam")
 
    # Step 1b: Mark duplicates on the coordinate-sorted BAM
    deduped_bam  = OUT_DIR / "marked_duplicates.bam"
    metrics_file = OUT_DIR / "duplicate_metrics.txt"
 
    run(gatk([
        "MarkDuplicates",
        "--INPUT",        sorted_bam,
        "--OUTPUT",       deduped_bam,
        "--METRICS_FILE", metrics_file,
        "--VALIDATION_STRINGENCY", "SILENT",
        "--CREATE_INDEX", "true",
    ]), "MarkDuplicates")
 
    print(f"\nStep 1 complete. Output: {deduped_bam}")
 
if __name__ == "__main__":
    main()
 