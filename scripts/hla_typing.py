"""
Step 3: HLA Typing

Two paths — clinical file takes priority, OptiType is the fallback:

  Path A — clinical HLA file (preferred):
    Patient already has HLA typed by targeted sequencing (Sanger SBT or
    NGS amplicon). Pass the file directly; no BAM processing needed.
    Accepts one allele per line, e.g.:

        HLA-A*31:01
        HLA-A*31:01
        HLA-B*37:01
        HLA-B*35:08
        HLA-C*06:02
        HLA-C*04:01

      python3 scripts/hla_typing.py --clinical path/to/hla.txt

  Path B — WES fallback (OptiType):
    No clinical file available. Extract HLA reads from the normal BAM
    (run create_test_data.sh first), then run OptiType.

    Phase 1 — remap + OptiType:
      python3 scripts/hla_typing.py --optitype

    Phase 2 — parse OptiType TSV → hla_alleles.txt:
      python3 scripts/hla_typing.py --parse

Output: results/step3/hla_alleles.txt
"""

import argparse
import csv
import re
import subprocess
import sys
from pathlib import Path

# ── Config ────────────────────────────────────────────────────────────────────

HLA_DIR       = Path.home() / "melanoma-pipeline/data/hla"
HLA_REF       = Path.home() / "melanoma-pipeline/reference/hla_reference_dna.fasta"
CANDIDATES_R1 = HLA_DIR / "hla_candidates_R1.fastq"
CANDIDATES_R2 = HLA_DIR / "hla_candidates_R2.fastq"
HLA_FISHED_R1 = HLA_DIR / "hla_fished_R1.fastq"
HLA_FISHED_R2 = HLA_DIR / "hla_fished_R2.fastq"
OUT_DIR       = Path.home() / "melanoma-pipeline/results/step3"
SAMPLE_PREFIX = "hcc1143_normal"

# ── Helpers ───────────────────────────────────────────────────────────────────

def run(cmd, step_name):
    print(f"\n[{step_name}] running...")
    result = subprocess.run([str(c) for c in cmd])
    if result.returncode != 0:
        print(f"FAILED: {step_name}")
        sys.exit(1)
    print(f"DONE: {step_name}")


def write_alleles(alleles, out_file, source):
    out_file.parent.mkdir(parents=True, exist_ok=True)
    with open(out_file, "w") as f:
        for a in alleles:
            f.write(a + "\n")
    print(f"\n── HLA Alleles ({source}) ──")
    for a in alleles:
        print(f"  {a}")
    print(f"\nStep 3 complete. Alleles saved to: {out_file}")

# ── Path A: clinical HLA file ─────────────────────────────────────────────────

# Accepted formats (all normalised to HLA-A*31:01):
#   HLA-A*31:01          (standard)
#   A*31:01              (no HLA- prefix)
#   A*31:01:02           (4-field — truncated to 2-field for pVACseq compat)
#   HLA-A*31:01:02:03

_ALLELE_RE = re.compile(
    r"(?:HLA-)?([ABC])\*(\d+):(\d+)",
    re.IGNORECASE,
)

def parse_clinical_file(path):
    """Parse a clinical HLA text file into normalised 2-field allele strings."""
    alleles = []
    with open(path) as f:
        for lineno, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            m = _ALLELE_RE.search(line)
            if not m:
                print(f"WARNING: line {lineno} not recognised as HLA allele: {line!r}")
                continue
            gene, field1, field2 = m.group(1).upper(), m.group(2), m.group(3)
            alleles.append(f"HLA-{gene}*{field1}:{field2}")

    if not alleles:
        print(f"ERROR: no valid HLA alleles found in {path}")
        sys.exit(1)

    return alleles

# ── Path B: OptiType fallback ─────────────────────────────────────────────────

def ensure_hla_ref_indexed(hla_ref):
    if not Path(str(hla_ref) + ".bwt").exists():
        run(["bwa", "index", str(hla_ref)], "Index HLA reference")
    else:
        print("[Index HLA reference] already indexed, skipping.")


def remap_to_hla_ref(hla_ref, r1, r2, out_r1, out_r2):
    """Remap candidate reads against HLA reference, split mapped reads into R1/R2."""
    for f in [r1, r2]:
        if not f.exists():
            print(f"ERROR: {f} not found. Run create_test_data.sh first.")
            sys.exit(1)

    print("\n[Remap to HLA reference] running...")
    # bwa mem → name-sort → samtools fastq -1/-2 (3-stage pipeline)
    # Name-sort ensures samtools fastq can correctly assign flag 0x40/0x80
    # to R1/R2 when pairs are broken by -F 4 filtering.
    bwa = subprocess.Popen(
        ["bwa", "mem", "-t", "8", str(hla_ref), str(r1), str(r2)],
        stdout=subprocess.PIPE,
    )
    sort = subprocess.Popen(
        ["samtools", "sort", "-n", "-@", "8", "-"],
        stdin=bwa.stdout,
        stdout=subprocess.PIPE,
    )
    bwa.stdout.close()
    fastq = subprocess.Popen(
        [
            "samtools", "fastq", "-F", "4",
            "-1", str(out_r1),
            "-2", str(out_r2),
            "-0", "/dev/null",
            "-s", "/dev/null",
            "-",
        ],
        stdin=sort.stdout,
    )
    sort.stdout.close()
    fastq.communicate()
    sort.wait()
    bwa.wait()
    if bwa.returncode != 0 or sort.returncode != 0 or fastq.returncode != 0:
        print("FAILED: Remap to HLA reference")
        sys.exit(1)
    reads = sum(1 for _ in open(out_r1)) // 4
    print(f"DONE: Remap to HLA reference — {reads} read pairs fished")


def run_optitype(hla_r1, hla_r2, out_dir, sample_prefix):
    """Run OptiType via Docker against fished HLA reads (paired-end R1 + R2)."""
    run([
        "docker", "run", "--rm",
        "-v", f"{hla_r1.parent}:/input",
        "-v", f"{out_dir}:/output",
        "fred2/optitype",
        "-i", f"/input/{hla_r1.name}", f"/input/{hla_r2.name}",
        "--dna",
        "-o", "/output",
        "-p", sample_prefix,
        "-v",
    ], "OptiType")


def parse_optitype(out_dir, sample_prefix):
    """Parse OptiType TSV output into a list of HLA alleles."""
    result_file = out_dir / f"{sample_prefix}_result.tsv"
    if not result_file.exists():
        print(f"ERROR: {result_file} not found. Run --optitype first.")
        sys.exit(1)

    with open(result_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        row = next(reader)

    alleles = []
    for gene in ["A1", "A2", "B1", "B2", "C1", "C2"]:
        allele = row[gene]
        if allele:
            alleles.append(f"HLA-{allele}")
    return alleles

# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Step 3: HLA Typing — clinical file (preferred) or OptiType (fallback)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--clinical", metavar="FILE",
        help="Path A: use clinical HLA file (Sanger SBT / NGS amplicon) — skips BAM processing",
    )
    group.add_argument(
        "--optitype", action="store_true",
        help="Path B phase 1: remap candidates → fish reads → run OptiType",
    )
    group.add_argument(
        "--parse", action="store_true",
        help="Path B phase 2: parse OptiType TSV → hla_alleles.txt",
    )
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    alleles_file = OUT_DIR / "hla_alleles.txt"

    if args.clinical:
        print("Step 3 — Path A: clinical HLA file")
        clinical_path = Path(args.clinical)
        if not clinical_path.exists():
            print(f"ERROR: {clinical_path} not found.")
            sys.exit(1)
        alleles = parse_clinical_file(clinical_path)
        write_alleles(alleles, alleles_file, source=f"clinical file: {clinical_path.name}")

    elif args.optitype:
        print("Step 3 — Path B, Phase 1: remap + OptiType")
        ensure_hla_ref_indexed(HLA_REF)
        remap_to_hla_ref(HLA_REF, CANDIDATES_R1, CANDIDATES_R2, HLA_FISHED_R1, HLA_FISHED_R2)
        run_optitype(HLA_FISHED_R1, HLA_FISHED_R2, OUT_DIR, SAMPLE_PREFIX)
        print(f"\nDone. Results written to {OUT_DIR}")
        print("Next: python3 scripts/hla_typing.py --parse")

    elif args.parse:
        print("Step 3 — Path B, Phase 2: parsing OptiType results")
        alleles = parse_optitype(OUT_DIR, SAMPLE_PREFIX)
        write_alleles(alleles, alleles_file, source="OptiType (WES fallback)")


if __name__ == "__main__":
    main()
